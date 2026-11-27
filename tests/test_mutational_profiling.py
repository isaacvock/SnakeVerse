from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import pyBigWig
import pysam
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib"))
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from call_snps import expand_target_genomic_alleles  # noqa: E402
from mutprof import (  # noqa: E402
    anomaly_calling_selection,
    editing_track_targets,
    normalize_mutprof_library_type,
    parse_configured_mut_types,
    render_anomaly_calling_args,
)
from samples import load_samples, required_columns_for_assay  # noqa: E402
from tracks import coverage_strand_arg  # noqa: E402
from validation import validate_resolved_config  # noqa: E402


def make_sample(
    sample_id: str,
    library_type: str = "treated",
    strandedness: str = "reverse",
    unit_id: str = "L001",
) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "unit_id": unit_id,
        "fastq_1": "R1.fastq.gz",
        "fastq_2": "R2.fastq.gz",
        "condition": "K562",
        "replicate": "1",
        "strandedness": strandedness,
        "library_type": library_type,
        "target": "TDP43",
        "pool_id": "",
        "input_sample_id": "",
        "inline_barcodes": "",
        "sra_id": "",
        "sra_layout": "",
    }


def make_mutprof_config(fasta_path: str = "genome.fa") -> dict:
    return {
        "schema_version": 2,
        "project": {"name": "test", "run_name": "mp"},
        "assay": "mutational_profiling",
        "samples": "samples.tsv",
        "results_dir": "results/mp",
        "preprocessing": {"input_stage": "raw"},
        "reference": {
            "name": "test_genome",
            "fasta": fasta_path,
            "gtf": "genes.gtf",
            "indexes": {"star": ""},
        },
        "tools": {
            "fastqc": {"tool": "fastqc", "params": {}},
            "fastp": {"tool": "fastp", "params": {}},
            "star": {
                "tool": "star",
                "params": {
                    "align": {
                        "outSAMattributes": "NH HI AS NM MD",
                        "quantMode": "GeneCounts TranscriptomeSAM",
                    }
                },
            },
            "samtools": {"tool": "samtools", "params": {}},
            "featurecounts": {"tool": "featurecounts", "params": {}},
            "deeptools": {"tool": "deeptools", "params": {}},
            "multiqc": {"tool": "multiqc", "params": {}},
        },
        "modules": {
            "fastq_qc": {"enabled": True},
            "trimming": {
                "enabled": True,
                "tool": "fastp",
                "adapters": {"mode": "auto"},
            },
            "alignment": {
                "enabled": True,
                "tool": "star",
                "params": {"align": {"outFilterMismatchNoverReadLmax": 0.08}},
            },
            "mark_duplicates": {"enabled": True},
            "bam_filter": {
                "enabled": True,
                "min_mapq": 10,
                "excluded_flags": 3852,
                "keep_duplicates": False,
            },
            "umi_extraction": {"enabled": False},
            "umi_deduplication": {"enabled": False},
            "bam_qc": {"enabled": True},
            "gene_counts": {
                "enabled": True,
                "tool": "featurecounts",
                "modes": {"standard": True},
                "strandedness": "sample",
            },
            "coverage": {"enabled": True, "strand_mode": "sample"},
            "anomaly_calling": {
                "enabled": True,
                "mode": "auto",
                "min_base_quality": 30,
                "min_depth": 5,
                "min_alt": 3,
                "min_freq": 0.25,
                "target_min_freq": 0.90,
                "err_rate": 0.002,
                "p_cutoff": 1e-4,
            },
            "mutation_counts": {
                "enabled": True,
                "mut_types": ["AG", "CT"],
                "min_base_quality": 30,
                "min_end_distance": 5,
                "indel_flank": 3,
                "min_depth": 1,
            },
            "editing_tracks": {
                "enabled": True,
                "min_depth": 5,
                "min_edit_count": 1,
                "emit_bedgraph": False,
            },
            "multiqc": {"enabled": True},
        },
    }


def test_mut_types_and_target_allele_expansion() -> None:
    cfg = {"modules": {"mutation_counts": {"mut_types": ["ag", "CT", "AG"]}}}
    assert parse_configured_mut_types(cfg) == ["AG", "CT"]
    assert expand_target_genomic_alleles("AG") == frozenset({("A", "G"), ("T", "C")})
    assert expand_target_genomic_alleles("AG,CT") == frozenset(
        {("A", "G"), ("T", "C"), ("C", "T"), ("G", "A")}
    )

    with pytest.raises(ValueError, match="Invalid mutation type"):
        parse_configured_mut_types({"modules": {"mutation_counts": {"mut_types": ["AA"]}}})
    with pytest.raises(ValueError, match="at least one"):
        parse_configured_mut_types({"modules": {"mutation_counts": {"mut_types": []}}})


def test_sample_sheet_loading_and_library_type_normalization(tmp_path: Path) -> None:
    assert "library_type" in required_columns_for_assay("mutational_profiling")
    assert normalize_mutprof_library_type("Negative Control") == "negative_control"

    tsv = tmp_path / "samples.tsv"
    tsv.write_text(
        "sample_id\tunit_id\tfastq_1\tcondition\treplicate\tstrandedness\tlibrary_type\n"
        "S1\tL1\ta.fq.gz\tWT\t1\tunstranded\ttreated\n"
        "C1\tL1\tb.fq.gz\tWT\t1\tunstranded\tnegative control\n",
        encoding="utf-8",
    )
    rows = load_samples(tsv, "mutational_profiling")
    assert rows[1]["library_type"] == "negative_control"
    assert (
        coverage_strand_arg(
            {"assay": "mutational_profiling", "modules": {"coverage": {"strand_mode": "sample"}}},
            [make_sample("S1", strandedness="reverse")],
            "S1",
        )
        == "--filterRNAstrand reverse"
    )


def test_anomaly_mode_selection_and_validation() -> None:
    cfg = make_mutprof_config()
    with_controls = [
        make_sample("RBP_1", "treated"),
        make_sample("CTRL_1", "negative_control"),
        make_sample("CTRL_2", "negative_control"),
    ]
    selected, uses_target = anomaly_calling_selection(cfg, with_controls)
    assert selected == ["CTRL_1", "CTRL_2"]
    assert uses_target is False
    rendered_ctrl = render_anomaly_calling_args(cfg, with_controls)
    assert "--target-mut-types" not in rendered_ctrl

    treated_only = [make_sample("RBP_1", "treated"), make_sample("RBP_2", "treated")]
    selected_t, uses_target_t = anomaly_calling_selection(cfg, treated_only)
    assert selected_t == ["RBP_1", "RBP_2"]
    assert uses_target_t is True
    rendered_treated = render_anomaly_calling_args(cfg, treated_only)
    assert "--target-mut-types AG,CT" in rendered_treated
    assert "--target-min-freq 0.9" in rendered_treated

    errors, _ = validate_resolved_config(cfg, with_controls)
    assert errors == []

    # Controls explicitly requested when none exist -> error
    bad_ctrl_cfg = deepcopy(cfg)
    bad_ctrl_cfg["modules"]["anomaly_calling"]["mode"] = "controls"
    errs, _ = validate_resolved_config(bad_ctrl_cfg, treated_only)
    assert any("no negative_control samples" in e for e in errs)

    # STAR missing MD tag -> error
    no_md_cfg = deepcopy(cfg)
    no_md_cfg["tools"]["star"]["params"]["align"]["outSAMattributes"] = "NH HI AS NM"
    errs_md, _ = validate_resolved_config(no_md_cfg, with_controls)
    assert any("MD" in e for e in errs_md)

    # editing_tracks without mutation_counts -> error
    no_mc_cfg = deepcopy(cfg)
    no_mc_cfg["modules"]["mutation_counts"]["enabled"] = False
    errs_et, _ = validate_resolved_config(no_mc_cfg, with_controls)
    assert any("editing_tracks requires modules.mutation_counts" in e for e in errs_et)

    # Disabling editing_tracks removes track targets cleanly
    tracks_on = editing_track_targets(cfg, [make_sample("RBP_1")], "results/mp")
    assert len(tracks_on) == 8  # 2 mut_types * 2 metrics * 2 strands


def _write_indexed_bam(
    bam_path: Path,
    chrom: str,
    chrom_len: int,
    records: list[pysam.AlignedSegment],
    sort_by_qname: bool = False,
) -> None:
    unsorted = bam_path.with_suffix(".unsorted.bam")
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": chrom_len}]}
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out_sf:
        for rec in records:
            out_sf.write(rec)
    if sort_by_qname:
        subprocess.run(
            ["samtools", "sort", "-n", "-o", str(bam_path), str(unsorted)],
            check=True,
        )
    else:
        subprocess.run(
            ["samtools", "sort", "-o", str(bam_path), str(unsorted)],
            check=True,
        )
        subprocess.run(["samtools", "index", str(bam_path)], check=True)
    unsorted.unlink()


def _make_seg(
    qname: str,
    seq: str,
    start_0: int,
    md: str,
    flag: int = 0,
    cigar: str | None = None,
    quals: list[int] | None = None,
) -> pysam.AlignedSegment:
    seg = pysam.AlignedSegment()
    seg.query_name = qname
    seg.query_sequence = seq
    seg.flag = flag
    seg.reference_id = 0
    seg.reference_start = start_0
    seg.mapping_quality = 255
    seg.cigarstring = cigar or f"{len(seq)}M"
    seg.query_qualities = pysam.qualitystring_to_array(
        "".join(chr(q + 33) for q in (quals or [37] * len(seq)))
    )
    seg.set_tag("MD", md, value_type="Z")
    return seg


def test_call_snps_pooled_controls_vs_mutation_type_aware_treated(tmp_path: Path) -> None:
    ref_seq = "A" * 100 + "T" * 100 + "C" * 100
    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">chrTest\n{ref_seq}\n", encoding="ascii")
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)

    # 1) Pooled negative controls: each control has 2 G's out of 4 reads at pos 21 (0-based 20).
    # Individually alt=2 < min_alt=3, but pooled depth=8 >= 5, alt=4 >= 3, freq=0.5 >= 0.25.
    for idx in (1, 2):
        recs = []
        for r_i in range(4):
            seq = list("A" * 40)
            md = "40"
            if r_i < 2:
                seq[20] = "G"
                md = "20A19"
            recs.append(_make_seg(f"c{idx}_{r_i}", "".join(seq), 0, md))
        _write_indexed_bam(tmp_path / f"ctrl{idx}.bam", "chrTest", len(ref_seq), recs)

    out_ctrl_txt = tmp_path / "ctrl_snps.txt"
    out_ctrl_tsv = tmp_path / "ctrl_snps.tsv"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "workflow" / "scripts" / "call_snps.py"),
            "-f",
            str(fasta),
            "-o",
            str(out_ctrl_txt),
            "--output-tsv",
            str(out_ctrl_tsv),
            "--min-depth",
            "5",
            "--min-alt",
            "3",
            "--min-freq",
            "0.25",
            "-b",
            str(tmp_path / "ctrl1.bam"),
            str(tmp_path / "ctrl2.bam"),
        ],
        check=True,
    )
    assert out_ctrl_txt.read_text(encoding="utf-8").strip().splitlines() == [
        "A:G:chrTest:21"
    ]

    # 2) Treated BAM with mutation-type-aware thresholds (AG -> A->G and T->C @ 0.90):
    # - pos 21 (0-based 20, ref A): 5/10 G (50% true + strand AG edit -> MUST NOT be masked)
    # - pos 31 (0-based 30, ref A): 10/10 G (100% homozygous SNP -> MUST be masked)
    # - pos 121 (0-based 120, ref T): 5/10 C (50% true - strand AG edit -> MUST NOT be masked)
    # - pos 221 (0-based 220, ref C): 5/10 A (50% non-target C->A SNP -> MUST be masked)
    t_recs = []
    for r_i in range(10):
        s1 = list("A" * 50)
        s1[30] = "G"
        if r_i < 5:
            s1[20] = "G"
            md1 = "20A9A19"
        else:
            md1 = "30A19"
        t_recs.append(_make_seg(f"t_a_{r_i}", "".join(s1), 0, md1))

        s2 = list("T" * 50)
        md2 = "50"
        if r_i < 5:
            s2[20] = "C"
            md2 = "20T29"
        t_recs.append(_make_seg(f"t_t_{r_i}", "".join(s2), 100, md2))

        s3 = list("C" * 50)
        md3 = "50"
        if r_i < 5:
            s3[20] = "A"
            md3 = "20C29"
        t_recs.append(_make_seg(f"t_c_{r_i}", "".join(s3), 200, md3))

    treated_bam = tmp_path / "treated.bam"
    _write_indexed_bam(treated_bam, "chrTest", len(ref_seq), t_recs)

    out_t_txt = tmp_path / "treated_snps.txt"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "workflow" / "scripts" / "call_snps.py"),
            "-f",
            str(fasta),
            "-o",
            str(out_t_txt),
            "--target-mut-types",
            "AG",
            "--target-min-freq",
            "0.90",
            "--min-freq",
            "0.25",
            "-b",
            str(treated_bam),
        ],
        check=True,
    )
    called = set(out_t_txt.read_text(encoding="utf-8").strip().splitlines())
    assert called == {"A:G:chrTest:31", "C:A:chrTest:221"}


def test_count_mutations_and_editing_tracks_stranded_and_unstranded(tmp_path: Path) -> None:
    # Reference: 0..99 = 'A', 100..199 = 'T' (so minus-strand RNA 100..199 is 'A')
    ref_seq = "A" * 100 + "T" * 100
    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">chrTest\n{ref_seq}\n", encoding="ascii")
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)

    # Create reverse-stranded PE reads:
    # - 6 overlapping pairs on '+' transcript (R1 reverse flag 99/83 -> R1 antisense: flag 0x1|0x2|0x40|0x10 = 83,
    #   R2 sense: flag 0x1|0x2|0x80 = 131... wait: proper reverse-stranded + transcript:
    #   R1 is reverse (0x10, flag=83), R2 is forward (flag=163)):
    #   Both mates cover 0-based 25 (1-based pos 26) with A->G in 3 of 6 fragments -> depth MUST be 6 (not 12!), G=3, A=3.
    #   Also put a SNP at 1-based pos 27 (0-based 26) masked via snps.txt, and a read-end mismatch at offset 1 (0-based 11 -> masked by min_end_distance=5).
    # - 6 single/paired reads on '-' transcript (100..160): R1 forward (flag=99), R2 reverse (flag=147).
    #   At 0-based 125 (1-based pos 126), + ref is 'T' -> '-' RNA ref is 'A'; BAM base 'C' -> '-' RNA base 'G' in 2 of 6 fragments!
    recs = []
    for i in range(6):
        s_r2 = list("A" * 40)  # covers 10..49
        s_r1 = list("A" * 40)  # covers 15..54 (overlaps 15..49, including 25 & 26!)
        # Read-end mismatch on R2 at qpos=1 (0-based 11, end_dist=1 < 5 -> masked)
        s_r2[1] = "G"
        # Masked SNP at 0-based 26 (pos 27)
        s_r2[16] = "C"
        s_r1[11] = "C"
        if i < 3:
            s_r2[15] = "G"  # 0-based 25 (pos 26)
            s_r1[10] = "G"  # 0-based 25 (pos 26)
            md_r2 = "1A13A0A23"
            md_r1 = "10A0A28"
        else:
            md_r2 = "1A14A23"
            md_r1 = "11A28"
        recs.append(_make_seg(f"plus_frag_{i}", "".join(s_r1), 15, md_r1, flag=83))
        recs.append(_make_seg(f"plus_frag_{i}", "".join(s_r2), 10, md_r2, flag=163))

    # Duplicate fragment (0x400) that must be ignored even if present
    dup_seq = list("A" * 40)
    dup_seq[15] = "G"
    recs.append(_make_seg("dup_frag", "".join(dup_seq), 10, "15A24", flag=163 | 0x400))

    # Minus-strand transcript fragments covering 110..154
    for i in range(6):
        m_r1 = list("T" * 40)  # covers 110..149, forward on + ref -> minus RNA
        m_r2 = list("T" * 40)  # covers 115..154, reverse on + ref -> minus RNA
        if i < 2:
            m_r1[15] = "C"  # 0-based 125 (pos 126): T->C on + ref == A->G on - RNA!
            m_r2[10] = "C"
            md_m1 = "15T24"
            md_m2 = "10T29"
        else:
            md_m1 = "40"
            md_m2 = "40"
        recs.append(_make_seg(f"minus_frag_{i}", "".join(m_r1), 110, md_m1, flag=99))
        recs.append(_make_seg(f"minus_frag_{i}", "".join(m_r2), 115, md_m2, flag=147))

    qname_bam = tmp_path / "sample.qname.bam"
    _write_indexed_bam(qname_bam, "chrTest", len(ref_seq), recs, sort_by_qname=True)

    snp_file = tmp_path / "snps.txt"
    snp_file.write_text("A:C:chrTest:27\n", encoding="utf-8")

    perbase_stranded = tmp_path / "stranded.tsv.gz"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "workflow" / "scripts" / "count_mutations.py"),
            "-b",
            str(qname_bam),
            "-o",
            str(perbase_stranded),
            "--strandedness",
            "reverse",
            "--snps",
            str(snp_file),
            "--min-base-quality",
            "30",
            "--min-end-distance",
            "5",
            "-t",
            "2",
        ],
        check=True,
    )

    with gzip.open(perbase_stranded, "rt", encoding="utf-8") as fh:
        rows = {(int(r["pos"]), r["strand"]): r for r in csv.DictReader(fh, delimiter="\t")}

    # Pos 27 was in snps.txt -> position-wide masked!
    assert (27, "+") not in rows
    # Pos 12 (0-based 11) was at qpos=1 (< min_end_distance=5) and not covered inward -> omitted!
    assert (12, "+") not in rows
    # Pos 26 on '+' strand: 6 fragments (not 12!), 3 A and 3 G
    r26 = rows[(26, "+")]
    assert r26["ref_base"] == "A"
    assert int(r26["depth"]) == 6
    assert int(r26["A"]) == 3
    assert int(r26["G"]) == 3
    # Pos 126 on '-' strand: complemented to RNA sense -> ref_base='A', depth=6, A=4, G=2
    r126 = rows[(126, "-")]
    assert r126["ref_base"] == "A"
    assert int(r126["depth"]) == 6
    assert int(r126["A"]) == 4
    assert int(r126["G"]) == 2

    # Build stranded BigWig + bedGraph editing tracks
    track_dir = tmp_path / "tracks_stranded"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "workflow" / "scripts" / "build_editing_tracks.py"),
            "--perbase",
            str(perbase_stranded),
            "--fai",
            f"{fasta}.fai",
            "--strandedness",
            "reverse",
            "--mut-types",
            "AG",
            "--out-dir",
            str(track_dir),
            "--sample",
            "S1",
            "--min-depth",
            "5",
            "--min-edit-count",
            "1",
            "--emit-bedgraph",
        ],
        check=True,
    )
    with pyBigWig.open(str(track_dir / "S1.AG.edit_rate.plus.bw")) as bw:
        assert bw.values("chrTest", 25, 26)[0] == pytest.approx(0.5)
    with pyBigWig.open(str(track_dir / "S1.AG.edit_count.minus.bw")) as bw:
        assert bw.values("chrTest", 125, 126)[0] == pytest.approx(2.0)
    assert (track_dir / "S1.AG.edit_rate.plus.bedGraph.gz").exists()

    # Unstranded mode on the exact same BAM emits +-reference strand counts with strand='.'
    perbase_unstranded = tmp_path / "unstranded.tsv.gz"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "workflow" / "scripts" / "count_mutations.py"),
            "-b",
            str(qname_bam),
            "-o",
            str(perbase_unstranded),
            "--strandedness",
            "unstranded",
            "--snps",
            str(snp_file),
        ],
        check=True,
    )
    with gzip.open(perbase_unstranded, "rt", encoding="utf-8") as fh:
        u_rows = {(int(r["pos"]), r["strand"]): r for r in csv.DictReader(fh, delimiter="\t")}
    # Pos 126 in unstranded mode stays on +-reference: ref_base='T', T=4, C=2, strand='.'
    u126 = u_rows[(126, ".")]
    assert u126["ref_base"] == "T"
    assert int(u126["T"]) == 4
    assert int(u126["C"]) == 2
