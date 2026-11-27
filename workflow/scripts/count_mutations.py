#!/usr/bin/env python3
"""Parallel per-nucleotide base-count and mutation caller for SnakeVerse.

Combines fastq2EZbakR's queryname-sorted BGZF parallel streaming, paired-end
mate-overlap reconciliation, indel-flank masking, read-end distance filtering,
and position-wide anomaly masking (mut_call.py) with a SnakeSHAPE/perbase-style
wide genomic nucleotide count output:

    chrom  pos  strand  ref_base  depth  A  C  G  T

Stranded libraries ('forward' or 'reverse') orient bases into the RNA sense
strand ('+' or '-'). Unstranded libraries ('unstranded') emit +-reference
strand base counts with strand='.'.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import heapq
import multiprocessing as mp
import os
import shutil
import struct
import sys
import tempfile
import zlib
from collections import defaultdict

import pysam


VALID_BASES = frozenset(("A", "C", "G", "T"))
BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
DNA_COMP = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "n": "n",
}
COMP_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")

SITE_SPILL_LIMIT = 500_000
MIN_PARALLEL_BAM_BYTES = 2 * 1024 * 1024

PERBASE_HEADER = ["chrom", "pos", "strand", "ref_base", "depth", "A", "C", "G", "T"]


def normalize_strandedness(raw: str) -> str:
    value = str(raw).strip().lower()
    mapping = {
        "forward": "forward",
        "f": "forward",
        "sense": "forward",
        "reverse": "reverse",
        "r": "reverse",
        "antisense": "reverse",
        "unstranded": "unstranded",
        "none": "unstranded",
        "u": "unstranded",
    }
    if value not in mapping:
        raise ValueError(
            f"Unsupported strandedness '{raw}'. Expected forward, reverse, or unstranded."
        )
    return mapping[value]


def load_snps(snp_path: str | None) -> dict[str, set[int]]:
    """Load REF:ALT:CHR:POS anomalies into chrom -> set of 0-based masked positions."""
    snp_by_chrom: dict[str, set[int]] = {}
    if not snp_path or not os.path.exists(snp_path):
        return snp_by_chrom

    with open(snp_path, "r", encoding="utf-8") as snp_file:
        for raw_line in snp_file:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split(":")
            if len(parts) < 4:
                continue
            chrom = ":".join(parts[2:-1])
            pos_0 = int(parts[-1]) - 1
            chrom_set = snp_by_chrom.get(chrom)
            if chrom_set is None:
                chrom_set = set()
                snp_by_chrom[chrom] = chrom_set
            chrom_set.add(pos_0)

    return snp_by_chrom


def _parse_bgzf_block_at(data: bytes, offset: int) -> tuple[int, bytes] | None:
    """Return `(block_len, uncompressed_bytes)` if `data[offset:]` starts a BGZF block."""
    n = len(data)
    if offset + 28 > n:
        return None
    if (
        data[offset] != 0x1F
        or data[offset + 1] != 0x8B
        or data[offset + 2] != 0x08
    ):
        return None
    flg = data[offset + 3]
    if (flg & 0x04) == 0 or (flg & 0xE0) != 0:
        return None

    xlen = struct.unpack_from("<H", data, offset + 10)[0]
    extra_end = offset + 12 + xlen
    if extra_end + 8 > n:
        return None

    bsize = None
    sub_pos = offset + 12
    while sub_pos + 4 <= extra_end:
        si1, si2, slen = struct.unpack_from("<BBH", data, sub_pos)
        if si1 == 66 and si2 == 67 and slen == 2 and sub_pos + 6 <= extra_end:
            bsize = struct.unpack_from("<H", data, sub_pos + 4)[0]
        sub_pos += 4 + slen

    if bsize is None:
        return None

    block_len = bsize + 1
    block_end = offset + block_len
    if block_end > n or block_len < 28:
        return None

    payload_start = extra_end
    if flg & 0x08:
        zero_idx = data.find(b"\x00", payload_start, block_end - 8)
        if zero_idx < 0:
            return None
        payload_start = zero_idx + 1
    if flg & 0x10:
        zero_idx = data.find(b"\x00", payload_start, block_end - 8)
        if zero_idx < 0:
            return None
        payload_start = zero_idx + 1
    if flg & 0x02:
        payload_start += 2

    if payload_start > block_end - 8:
        return None

    crc_expected, isize_expected = struct.unpack_from("<II", data, block_end - 8)
    if isize_expected > 65536:
        return None

    try:
        raw = zlib.decompress(data[payload_start : block_end - 8], -15)
    except zlib.error:
        return None

    if len(raw) != isize_expected:
        return None
    if (zlib.crc32(raw) & 0xFFFFFFFF) != crc_expected:
        return None

    return block_len, raw


def _is_valid_bam_record(ubuf: bytes, pos: int, ref_lengths: list[int]) -> int | str | None:
    """Validate whether `ubuf[pos:]` begins a binary BAM alignment record."""
    n = len(ubuf)
    if pos + 36 > n:
        return "TRUNCATED"

    (
        block_size,
        ref_id,
        aln_pos,
        l_read_name,
        _mapq,
        _bin,
        n_cigar_op,
        flag,
        l_seq,
        next_ref_id,
        next_pos,
        _tlen,
    ) = struct.unpack_from("<iiiBBHHHiiii", ubuf, pos)

    n_refs = len(ref_lengths)
    if l_read_name < 2 or l_seq < 0 or (flag & 0xF000) != 0:
        return None
    if ref_id < -1 or ref_id >= n_refs or next_ref_id < -1 or next_ref_id >= n_refs:
        return None
    if aln_pos < -1 or (ref_id >= 0 and aln_pos >= ref_lengths[ref_id]):
        return None
    if next_pos < -1 or (next_ref_id >= 0 and next_pos >= ref_lengths[next_ref_id]):
        return None

    min_body = 32 + l_read_name + 4 * n_cigar_op + ((l_seq + 1) // 2) + l_seq
    if block_size < min_body or block_size > 64 * 1024 * 1024:
        return None

    rec_end = pos + 4 + block_size
    cigar_end = pos + 36 + l_read_name + 4 * n_cigar_op
    if cigar_end > n:
        return "TRUNCATED"

    rname_bytes = ubuf[pos + 36 : pos + 36 + l_read_name]
    if rname_bytes[-1] != 0:
        return None
    for b in rname_bytes[:-1]:
        if b < 33 or b > 126:
            return None

    if n_cigar_op > 0:
        q_sum = 0
        cigar_offset = pos + 36 + l_read_name
        for idx in range(n_cigar_op):
            cigar_int = struct.unpack_from("<I", ubuf, cigar_offset + 4 * idx)[0]
            op = cigar_int & 0xF
            oplen = cigar_int >> 4
            if op > 8 or oplen <= 0:
                return None
            if op in (0, 1, 4, 7, 8):
                q_sum += oplen
        if l_seq > 0 and q_sum != l_seq:
            return None

    if rec_end > n:
        return "TRUNCATED"
    return rec_end


def _find_voffset_near(
    fh, file_size: int, min_block_offset: int, target_byte: int, ref_lengths: list[int]
) -> int | None:
    """Find a verified BAM record virtual offset `(block_offset << 16) | uoffset`."""
    search_pos = max(min_block_offset, target_byte)
    if search_pos >= file_size - 28:
        return None

    fh.seek(search_pos)
    window = fh.read(min(512 * 1024, file_size - search_pos))

    b0_rel = None
    scan = 0
    limit = len(window) - 28
    while scan <= limit:
        idx = window.find(b"\x1f\x8b\x08", scan)
        if idx < 0 or idx > limit:
            break
        parsed0 = _parse_bgzf_block_at(window, idx)
        if parsed0 is not None:
            blen0, _ = parsed0
            next1 = idx + blen0
            if search_pos + next1 == file_size:
                b0_rel = idx
                break
            parsed1 = _parse_bgzf_block_at(window, next1)
            if parsed1 is not None:
                blen1, _ = parsed1
                next2 = next1 + blen1
                if search_pos + next2 == file_size or _parse_bgzf_block_at(window, next2) is not None:
                    b0_rel = idx
                    break
        scan = idx + 1

    if b0_rel is None:
        return None

    cur_file_offset = search_pos + b0_rel
    while cur_file_offset < file_size - 28:
        fh.seek(cur_file_offset)
        buf = fh.read(min(512 * 1024, file_size - cur_file_offset))
        blocks = []
        rel = 0
        reached_eof = False
        while rel < len(buf) and len(blocks) < 6:
            parsed = _parse_bgzf_block_at(buf, rel)
            if parsed is None:
                break
            blen, raw = parsed
            if len(raw) == 0 and cur_file_offset + rel + blen == file_size:
                reached_eof = True
                break
            blocks.append((cur_file_offset + rel, blen, raw))
            rel += blen
            if cur_file_offset + rel == file_size:
                reached_eof = True
                break

        if not blocks:
            return None

        first_block_offset, first_blen, first_raw = blocks[0]
        if len(first_raw) == 0:
            cur_file_offset = first_block_offset + first_blen
            continue

        ubuf = b"".join(b[2] for b in blocks)
        first_u_len = len(first_raw)

        for uoff in range(first_u_len):
            cur_p = uoff
            valid_chain = 0
            ok = True
            while valid_chain < 4:
                nxt = _is_valid_bam_record(ubuf, cur_p, ref_lengths)
                if nxt is None:
                    ok = False
                    break
                if nxt == "TRUNCATED":
                    ok = valid_chain >= 2
                    break
                valid_chain += 1
                cur_p = nxt
                if cur_p == len(ubuf) and reached_eof:
                    break
            if ok and valid_chain >= 1:
                return (first_block_offset << 16) | uoff

        cur_file_offset = first_block_offset + first_blen

    return None


def compute_bam_chunks(
    bam_path: str, threads: int
) -> tuple[list[tuple[int, int | None]], list[str], bool]:
    """Return `(chunks, reference_names, is_empty)` across `bam_path`."""
    with pysam.AlignmentFile(bam_path, "rb") as sf:
        references = list(sf.references)
        first_voff = sf.tell()
        try:
            next(sf)
        except StopIteration:
            return [], references, True
        ref_lengths = list(sf.lengths)

    file_size = os.path.getsize(bam_path)
    if threads <= 1 or file_size < MIN_PARALLEL_BAM_BYTES:
        return [(first_voff, None)], references, False

    min_block_offset = (first_voff >> 16) + 1
    split_voffs = [first_voff]

    with open(bam_path, "rb") as fh, pysam.AlignmentFile(bam_path, "rb") as verify_sf:
        for k in range(1, threads):
            target = (k * file_size) // threads
            voff = _find_voffset_near(fh, file_size, min_block_offset, target, ref_lengths)
            if voff is None or voff <= split_voffs[-1]:
                continue
            try:
                verify_sf.seek(voff)
                rec = next(verify_sf)
                if not rec.query_name:
                    continue
            except (OSError, StopIteration, ValueError):
                continue
            split_voffs.append(voff)

    chunks = []
    for idx, start_voff in enumerate(split_voffs):
        end_voff = split_voffs[idx + 1] if idx + 1 < len(split_voffs) else None
        chunks.append((start_voff, end_voff))
    return chunks, references, False


def infer_record_strand(rec: pysam.AlignedSegment, strandedness: str) -> str:
    """Return '+', '-', or '.' for one aligned segment."""
    if strandedness == "unstranded":
        return "."
    strand_is_forward = strandedness == "forward"
    if rec.is_paired:
        if rec.is_read2 and not rec.is_read1:
            same_as_rna = not strand_is_forward
        else:
            same_as_rna = strand_is_forward
    else:
        same_as_rna = strand_is_forward
    aligned_strand = "-" if rec.is_reverse else "+"
    if same_as_rna:
        return aligned_strand
    return "+" if aligned_strand == "-" else "-"


def _flush_site_spill(site_counts: dict, spill_dir: str, worker_idx: int, spill_idx: int) -> str:
    """Sort `site_counts` by `(ref_id, pos_1based, strand)` and spill to TSV."""
    spill_path = os.path.join(spill_dir, f"sites_{worker_idx}_{spill_idx}.tsv")
    with open(spill_path, "w", newline="", encoding="utf-8") as fh:
        for (ref_id, pos_1, strand), (ref_base, a_c, c_c, g_c, t_c) in sorted(
            site_counts.items()
        ):
            fh.write(f"{ref_id}\t{pos_1}\t{strand}\t{ref_base}\t{a_c}\t{c_c}\t{g_c}\t{t_c}\n")
    site_counts.clear()
    return spill_path


def _iter_site_spill(spill_path: str):
    """Yield `((ref_id, pos_1, strand), ref_base, a, c, g, t)` from a sorted spill."""
    with open(spill_path, "r", encoding="utf-8") as fh:
        for line in fh:
            rid_s, pos_s, strand, ref_base, a_s, c_s, g_s, t_s = line.rstrip("\n").split("\t")
            yield (
                (int(rid_s), int(pos_s), strand),
                ref_base,
                int(a_s),
                int(c_s),
                int(g_s),
                int(t_s),
            )


def _merge_site_spills(
    spill_paths: list[str],
    references: list[str],
    output_path: str,
    min_depth: int,
) -> int:
    """Stream-merge sorted site spill files into the final wide per-base TSV[.gz]."""
    is_gz = output_path.endswith(".gz")
    opener = (
        (lambda p: gzip.open(p, "wt", newline="", encoding="utf-8", compresslevel=3))
        if is_gz
        else (lambda p: open(p, "w", newline="", encoding="utf-8"))
    )

    streams = [_iter_site_spill(p) for p in spill_paths]
    merged = heapq.merge(*streams, key=lambda item: item[0])

    emitted = 0
    with opener(output_path) as out_fh:
        wr = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        wr.writerow(PERBASE_HEADER)

        cur_key = None
        cur_ref = ""
        cur_a = 0
        cur_c = 0
        cur_g = 0
        cur_t = 0

        for key, ref_base, a_c, c_c, g_c, t_c in merged:
            if key == cur_key:
                cur_a += a_c
                cur_c += c_c
                cur_g += g_c
                cur_t += t_c
            else:
                if cur_key is not None:
                    depth = cur_a + cur_c + cur_g + cur_t
                    if depth >= min_depth:
                        wr.writerow(
                            [
                                references[cur_key[0]],
                                cur_key[1],
                                cur_key[2],
                                cur_ref,
                                depth,
                                cur_a,
                                cur_c,
                                cur_g,
                                cur_t,
                            ]
                        )
                        emitted += 1
                cur_key = key
                cur_ref = ref_base
                cur_a = a_c
                cur_c = c_c
                cur_g = g_c
                cur_t = t_c

        if cur_key is not None:
            depth = cur_a + cur_c + cur_g + cur_t
            if depth >= min_depth:
                wr.writerow(
                    [
                        references[cur_key[0]],
                        cur_key[1],
                        cur_key[2],
                        cur_ref,
                        depth,
                        cur_a,
                        cur_c,
                        cur_g,
                        cur_t,
                    ]
                )
                emitted += 1

    return emitted


def _process_subfragment(
    records: list[pysam.AlignedSegment],
    ref_id: int,
    rname: str,
    strand: str,
    min_base_quality: int,
    min_end_distance: int,
    indel_flank: int,
    chrom_snps: set[int] | None,
    site_counts: dict,
) -> None:
    """Reconcile overlapping mates for one `(qname, ref_id, strand)` group and accumulate counts."""
    complement = strand == "-"
    md_store: dict[int, tuple[bool, int, int, bool, str, str]] = {}
    indel_masked_rpos: set[int] = set()

    for r in records:
        cigar_str = r.cigarstring
        if cigar_str is None or r.query_sequence is None:
            continue

        bad_qpos: set[int] | None = None
        if indel_flank >= 0 and ("I" in cigar_str or "D" in cigar_str):
            bad_qpos = set()
            q_cur = 0
            r_cur = r.reference_start
            for op, length in r.cigartuples or ():
                if op in (0, 7, 8):  # M, =, X
                    q_cur += length
                    r_cur += length
                elif op == 1:  # I
                    for qp in range(max(0, q_cur - indel_flank), q_cur + length + indel_flank):
                        bad_qpos.add(qp)
                    for rp in range(r_cur - indel_flank, r_cur + indel_flank):
                        indel_masked_rpos.add(rp)
                    q_cur += length
                elif op == 2:  # D
                    for qp in range(max(0, q_cur - indel_flank), q_cur + indel_flank):
                        bad_qpos.add(qp)
                    for rp in range(r_cur - indel_flank, r_cur + length + indel_flank):
                        indel_masked_rpos.add(rp)
                    r_cur += length
                elif op == 3:  # N (intron splice skip)
                    r_cur += length
                elif op == 4:  # S (soft clip)
                    q_cur += length

        try:
            aligned_pairs = r.get_aligned_pairs(matches_only=True, with_seq=True)
        except ValueError as exc:
            raise RuntimeError(
                f"Alignment '{r.query_name}' on '{rname}' is missing a valid MD tag "
                f"required for mutation calling: {exc}"
            ) from exc

        seq = r.query_sequence
        if complement:
            seq = seq.translate(COMP_TABLE)
        quals = r.query_qualities
        if quals is None:
            quals = [60] * len(seq)

        q_start = r.query_alignment_start
        q_len = r.query_alignment_length

        for qpos, rpos, rbase in aligned_pairs:
            if bad_qpos is not None and qpos in bad_qpos:
                continue
            if rbase is None:
                continue

            ref_oriented = DNA_COMP.get(rbase, "N") if complement else rbase
            ref_upper = ref_oriented.upper()
            if ref_upper not in VALID_BASES:
                continue

            read_upper = seq[qpos].upper()
            qual = int(quals[qpos])
            rel_q = qpos - q_start
            left_d = rel_q
            right_d = (q_len - 1) - rel_q
            end_dist = left_d if left_d <= right_d else right_d

            passes = (
                (read_upper in VALID_BASES)
                and (qual >= min_base_quality)
                and (end_dist >= min_end_distance)
            )
            obs = (
                passes,
                qual,
                end_dist,
                read_upper == ref_upper,
                ref_upper,
                read_upper,
            )

            prev = md_store.get(rpos)
            if prev is None or obs > prev:
                md_store[rpos] = obs

    if indel_masked_rpos:
        for pos in indel_masked_rpos:
            md_store.pop(pos, None)

    if not md_store:
        return

    for rpos, (passes, _qual, _end_dist, _is_ref, ref_upper, read_upper) in md_store.items():
        if not passes:
            continue
        if chrom_snps is not None and rpos in chrom_snps:
            continue

        key = (ref_id, rpos + 1, strand)
        entry = site_counts.get(key)
        if entry is None:
            entry = [ref_upper, 0, 0, 0, 0]
            site_counts[key] = entry
        entry[1 + BASE_TO_IDX[read_upper]] += 1


def _process_read_group(
    records: list[pysam.AlignedSegment],
    strandedness: str,
    min_base_quality: int,
    min_end_distance: int,
    indel_flank: int,
    snp_by_chrom: dict[str, set[int]],
    site_counts: dict,
) -> None:
    """Filter primary records sharing one `query_name` and accumulate per-base counts."""
    grouped: dict[tuple[int, str, str], list[pysam.AlignedSegment]] = defaultdict(list)
    for r in records:
        if (
            r.is_unmapped
            or r.is_secondary
            or r.is_qcfail
            or r.is_duplicate
            or r.is_supplementary
            or r.reference_id < 0
            or r.reference_name is None
            or r.cigarstring is None
            or r.query_sequence is None
        ):
            continue
        strand = infer_record_strand(r, strandedness)
        grouped[(r.reference_id, r.reference_name, strand)].append(r)

    for (ref_id, rname, strand), sub_records in grouped.items():
        _process_subfragment(
            sub_records,
            ref_id,
            rname,
            strand,
            min_base_quality,
            min_end_distance,
            indel_flank,
            snp_by_chrom.get(rname),
            site_counts,
        )


_WORKER_CFG: dict = {}


def _init_worker(cfg: dict) -> None:
    global _WORKER_CFG
    _WORKER_CFG = cfg


def _run_chunk_worker(task: tuple[int, int, int | None]) -> list[str]:
    """Process BAM records in `[start_voff, end_voff]` with exact `query_name` handoff."""
    worker_idx, start_voff, end_voff = task
    cfg = _WORKER_CFG

    bam_path = cfg["bam_path"]
    strandedness = cfg["strandedness"]
    min_base_quality = cfg["min_base_quality"]
    min_end_distance = cfg["min_end_distance"]
    indel_flank = cfg["indel_flank"]
    snp_by_chrom = cfg["snp_by_chrom"]
    spill_dir = cfg["spill_dir"]

    site_counts: dict = {}
    spills: list[str] = []
    spill_idx = 0

    with pysam.AlignmentFile(bam_path, "rb") as sf:
        sf.seek(start_voff)

        first_rec = None
        first_rec_voff = None
        if worker_idx > 0:
            try:
                r0 = next(sf)
            except StopIteration:
                return spills
            skip_qname = r0.query_name
            while True:
                rec_voff = sf.tell()
                try:
                    r = next(sf)
                except StopIteration:
                    return spills
                if r.query_name != skip_qname:
                    if end_voff is not None and rec_voff > end_voff:
                        return spills
                    first_rec = r
                    first_rec_voff = rec_voff
                    break
        else:
            rec_voff = sf.tell()
            try:
                first_rec = next(sf)
                first_rec_voff = rec_voff
            except StopIteration:
                return spills

        cur_qname = first_rec.query_name
        cur_group = [first_rec]
        passed_end = end_voff is not None and first_rec_voff >= end_voff

        while True:
            rec_voff = sf.tell()
            try:
                r = next(sf)
            except StopIteration:
                break

            qname = r.query_name
            if qname == cur_qname:
                cur_group.append(r)
                if end_voff is not None and rec_voff >= end_voff:
                    passed_end = True
            else:
                _process_read_group(
                    cur_group,
                    strandedness,
                    min_base_quality,
                    min_end_distance,
                    indel_flank,
                    snp_by_chrom,
                    site_counts,
                )
                if len(site_counts) >= SITE_SPILL_LIMIT:
                    spills.append(_flush_site_spill(site_counts, spill_dir, worker_idx, spill_idx))
                    spill_idx += 1

                if passed_end or (end_voff is not None and rec_voff > end_voff):
                    cur_group = None
                    break

                cur_qname = qname
                cur_group = [r]
                if end_voff is not None and rec_voff >= end_voff:
                    passed_end = True

        if cur_group:
            _process_read_group(
                cur_group,
                strandedness,
                min_base_quality,
                min_end_distance,
                indel_flank,
                snp_by_chrom,
                site_counts,
            )

    if site_counts:
        spills.append(_flush_site_spill(site_counts, spill_dir, worker_idx, spill_idx))

    return spills


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Tabulate strand-aware genomic per-nucleotide base counts from a queryname-sorted BAM"
    )
    parser.add_argument(
        "-b",
        "--bam",
        type=str,
        required=True,
        help="Queryname-sorted BAM file to process (must contain MD tags)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output wide per-base TSV or TSV.GZ path",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of parallel BGZF worker processes (default: 1)",
    )
    parser.add_argument(
        "--strandedness",
        default="reverse",
        type=str,
        help="Library strandedness: forward, reverse, or unstranded (default: reverse)",
    )
    parser.add_argument(
        "--snps",
        default=None,
        help="Optional snps.txt file (REF:ALT:CHR:POS) for position-wide masking",
    )
    parser.add_argument(
        "--min-base-quality",
        default=30,
        type=int,
        help="Minimum Phred base quality required to count a base (default: 30)",
    )
    parser.add_argument(
        "--min-end-distance",
        default=5,
        type=int,
        help="Minimum distance in bp from either aligned read end (default: 5; 0 disables)",
    )
    parser.add_argument(
        "--indel-flank",
        default=3,
        type=int,
        help="Bases flanking any I/D CIGAR operation to mask (default: 3; <0 disables)",
    )
    parser.add_argument(
        "--min-depth",
        default=1,
        type=int,
        help="Minimum passing fragment depth required to emit a position (default: 1)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    strandedness = normalize_strandedness(args.strandedness)
    min_base_quality = max(0, int(args.min_base_quality))
    min_end_distance = max(0, int(args.min_end_distance))
    indel_flank = int(args.indel_flank)
    min_depth = max(1, int(args.min_depth))
    threads = max(1, int(args.threads))

    out_abs = os.path.abspath(args.output)
    out_dir = os.path.dirname(out_abs)
    os.makedirs(out_dir, exist_ok=True)

    print(
        f"Start count_mutations: {datetime.datetime.now()} "
        f"[strandedness={strandedness}, Phred>={min_base_quality}, "
        f"end_dist>={min_end_distance}bp, indel_flank={indel_flank}bp, min_depth>={min_depth}]",
        flush=True,
    )

    snp_by_chrom = load_snps(args.snps)
    masked_total = sum(len(v) for v in snp_by_chrom.values())
    print(f"Loaded {masked_total} masked genomic positions from {args.snps or '<none>'}.", flush=True)

    chunks, references, is_empty = compute_bam_chunks(args.bam, threads)
    if is_empty:
        _merge_site_spills([], references, out_abs, min_depth)
        print(f"Input BAM contains zero alignments; wrote header-only {out_abs}.", flush=True)
        return

    tmp_dir = tempfile.mkdtemp(prefix="count_mut_", dir=out_dir)
    try:
        cfg = {
            "bam_path": args.bam,
            "strandedness": strandedness,
            "min_base_quality": min_base_quality,
            "min_end_distance": min_end_distance,
            "indel_flank": indel_flank,
            "snp_by_chrom": snp_by_chrom,
            "spill_dir": tmp_dir,
        }

        tasks = [(idx, s_voff, e_voff) for idx, (s_voff, e_voff) in enumerate(chunks)]
        all_spills: list[str] = []
        if len(tasks) == 1:
            _init_worker(cfg)
            all_spills = _run_chunk_worker(tasks[0])
        else:
            with mp.Pool(processes=len(tasks), initializer=_init_worker, initargs=(cfg,)) as pool:
                for worker_spills in pool.map(_run_chunk_worker, tasks):
                    all_spills.extend(worker_spills)

        emitted = _merge_site_spills(all_spills, references, out_abs, min_depth)
        print(
            f"Finished count_mutations: {emitted} genomic site rows written to {out_abs} "
            f"({datetime.datetime.now()}).",
            flush=True,
        )
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR in count_mutations.py: {exc}", file=sys.stderr, flush=True)
        raise
