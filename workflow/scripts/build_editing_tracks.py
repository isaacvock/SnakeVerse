#!/usr/bin/env python3
"""Build stranded BigWig (and optional bedGraph.gz) tracks of editing rates and
edited read counts from a SnakeVerse per-base wide TSV[.gz] file.

For stranded samples ('forward' or 'reverse'), rows with strand=='+' populate
'.plus.' tracks and rows with strand=='-' populate '.minus.' tracks (since minus
rows in per-base TSV are already complemented into RNA sense orientation).

For unstranded samples ('unstranded', where rows have strand=='.' on the
+-reference strand), '.plus.' tracks report X->Y sites (ref_base == X) and
'.minus.' tracks report complementary comp(X)->comp(Y) sites corresponding to
minus-strand X->Y edits.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict

import pyBigWig


VALID_BASES = frozenset(("A", "C", "G", "T"))
DNA_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
BASE_COL_IDX = {"A": 5, "C": 6, "G": 7, "T": 8}


def read_fai_header(fai_path: str) -> list[tuple[str, int]]:
    header: list[tuple[str, int]] = []
    with open(fai_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2 and int(parts[1]) > 0:
                header.append((parts[0], int(parts[1])))
    if not header:
        raise ValueError(f"No valid contigs found in FASTA .fai: {fai_path}")
    return header


def parse_mut_types(raw: str) -> list[str]:
    result: list[str] = []
    for token in raw.split(","):
        mt = token.strip().upper()
        if not mt:
            continue
        if len(mt) != 2 or mt[0] not in VALID_BASES or mt[1] not in VALID_BASES or mt[0] == mt[1]:
            raise ValueError(
                f"Invalid mutation type '{token}'; expected two distinct ACGT bases such as AG or CT"
            )
        if mt not in result:
            result.append(mt)
    if not result:
        raise ValueError("At least one mutation type is required for editing tracks")
    return result


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def write_bigwig(
    path: str,
    chrom_header: list[tuple[str, int]],
    entries_by_chrom: dict[str, list[tuple[int, float]]],
) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    bw = pyBigWig.open(path, "w")
    try:
        bw.addHeader(chrom_header)
        for chrom, chrom_len in chrom_header:
            pts = entries_by_chrom.get(chrom)
            if not pts:
                continue
            pts.sort(key=lambda item: item[0])
            starts = [p[0] for p in pts if 0 <= p[0] < chrom_len]
            if not starts:
                continue
            ends = [s + 1 for s in starts]
            vals = [float(p[1]) for p in pts if 0 <= p[0] < chrom_len]
            bw.addEntries([chrom] * len(starts), starts, ends=ends, values=vals)
    finally:
        bw.close()


def write_bedgraph_gz(
    path: str,
    chrom_header: list[tuple[str, int]],
    entries_by_chrom: dict[str, list[tuple[int, float]]],
    is_integer: bool,
) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="\n", compresslevel=3) as fh:
        for chrom, chrom_len in chrom_header:
            pts = entries_by_chrom.get(chrom)
            if not pts:
                continue
            pts.sort(key=lambda item: item[0])
            for start, val in pts:
                if 0 <= start < chrom_len:
                    val_str = str(int(round(val))) if is_integer else f"{val:.6f}"
                    fh.write(f"{chrom}\t{start}\t{start + 1}\t{val_str}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate stranded BigWig (and optional bedGraph.gz) editing rate & count tracks"
    )
    parser.add_argument("--perbase", required=True, help="Input per-base wide TSV[.gz]")
    parser.add_argument("--fai", required=True, help="Reference FASTA .fai for contig lengths")
    parser.add_argument("--strandedness", required=True, help="Sample strandedness")
    parser.add_argument("--mut-types", required=True, help="Comma-separated mutation types, e.g. AG,CT")
    parser.add_argument("--out-dir", required=True, help="Output directory for editing tracks")
    parser.add_argument("--sample", required=True, help="Sample identifier prefix")
    parser.add_argument("--min-depth", type=int, default=5, help="Minimum depth to emit a site (default: 5)")
    parser.add_argument(
        "--min-edit-count",
        type=int,
        default=1,
        help="Minimum edited read count to emit a site (default: 1)",
    )
    parser.add_argument(
        "--emit-bedgraph",
        action="store_true",
        help="Also emit gzipped bedGraph (.bedGraph.gz) files alongside BigWigs",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    chrom_header = read_fai_header(args.fai)
    known_chroms = {c for c, _ in chrom_header}
    mut_types = parse_mut_types(args.mut_types)
    unstranded = args.strandedness.strip().lower() in {"unstranded", "none", "u"}
    min_depth = max(1, int(args.min_depth))
    min_edit_count = max(0, int(args.min_edit_count))

    # Map (mt, strand_label) -> {chrom: [(start_0, edit_rate), ...]} and edit_count
    rate_tracks: dict[tuple[str, str], dict[str, list[tuple[int, float]]]] = {
        (mt, sl): defaultdict(list) for mt in mut_types for sl in ("plus", "minus")
    }
    count_tracks: dict[tuple[str, str], dict[str, list[tuple[int, float]]]] = {
        (mt, sl): defaultdict(list) for mt in mut_types for sl in ("plus", "minus")
    }

    # Precompute lookup rules: (row_strand, ref_base) -> list of (mt, strand_label, alt_col_idx)
    lookup: dict[tuple[str, str], list[tuple[str, str, int]]] = defaultdict(list)
    if unstranded:
        for mt in mut_types:
            ref_b, alt_b = mt[0], mt[1]
            lookup[(".", ref_b)].append((mt, "plus", BASE_COL_IDX[alt_b]))
            lookup[(".", DNA_COMP[ref_b])].append(
                (mt, "minus", BASE_COL_IDX[DNA_COMP[alt_b]])
            )
    else:
        for mt in mut_types:
            ref_b, alt_b = mt[0], mt[1]
            lookup[("+", ref_b)].append((mt, "plus", BASE_COL_IDX[alt_b]))
            lookup[("-", ref_b)].append((mt, "minus", BASE_COL_IDX[alt_b]))

    with open_text(args.perbase) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError(f"Empty per-base file: {args.perbase}")
        for row in reader:
            if len(row) < 9:
                continue
            chrom = row[0]
            if chrom not in known_chroms:
                continue
            targets = lookup.get((row[2], row[3]))
            if not targets:
                continue
            depth = int(row[4])
            if depth < min_depth:
                continue
            start_0 = int(row[1]) - 1
            for mt, strand_label, alt_idx in targets:
                edit_cnt = int(row[alt_idx])
                if edit_cnt < min_edit_count:
                    continue
                rate_tracks[(mt, strand_label)][chrom].append((start_0, edit_cnt / depth))
                count_tracks[(mt, strand_label)][chrom].append((start_0, float(edit_cnt)))

    for mt in mut_types:
        for strand_label in ("plus", "minus"):
            rate_bw = os.path.join(
                args.out_dir, f"{args.sample}.{mt}.edit_rate.{strand_label}.bw"
            )
            count_bw = os.path.join(
                args.out_dir, f"{args.sample}.{mt}.edit_count.{strand_label}.bw"
            )
            write_bigwig(rate_bw, chrom_header, rate_tracks[(mt, strand_label)])
            write_bigwig(count_bw, chrom_header, count_tracks[(mt, strand_label)])

            if args.emit_bedgraph:
                rate_bg = os.path.join(
                    args.out_dir, f"{args.sample}.{mt}.edit_rate.{strand_label}.bedGraph.gz"
                )
                count_bg = os.path.join(
                    args.out_dir, f"{args.sample}.{mt}.edit_count.{strand_label}.bedGraph.gz"
                )
                write_bedgraph_gz(
                    rate_bg, chrom_header, rate_tracks[(mt, strand_label)], is_integer=False
                )
                write_bedgraph_gz(
                    count_bg, chrom_header, count_tracks[(mt, strand_label)], is_integer=True
                )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR in build_editing_tracks.py: {exc}", file=sys.stderr, flush=True)
        raise
