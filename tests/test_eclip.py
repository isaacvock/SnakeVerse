from __future__ import annotations

import gzip
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import pysam
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib"))

from clip import (  # noqa: E402
    crosslink_offset,
    dedup_args,
    eclip_barcode_adapter_overrides,
    representative_count_samtools_args,
    strand_flag_args,
    validate_qname_regex,
)
from params import render_tool_params, tool_extra, trimming_adapter_params  # noqa: E402
from samples import (  # noqa: E402
    alignment_fastq_1,
    alignment_fastq_2,
    bowtie2_reads_arg,
    bwa_mem2_reads_arg,
    enrichment_pairs,
    parse_inline_barcodes,
    processed_sample_layout,
    processed_unit_layout,
    star_reads_arg,
)
from validation import _validate_adapters, _validate_eclip, _validate_tool_extras  # noqa: E402


def sample(
    sample_id: str,
    library_type: str = "ip",
    control: str = "",
    unit_id: str = "L001",
    inline_barcodes: str = "",
) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "unit_id": unit_id,
        "fastq_1": "R1.fastq.gz",
        "fastq_2": "R2.fastq.gz",
        "condition": "K562",
        "replicate": "1",
        "library_type": library_type,
        "target": "PUM2",
        "pool_id": "",
        "input_sample_id": control,
        "inline_barcodes": inline_barcodes,
        "sra_id": "",
        "sra_layout": "",
    }


def eclip_config() -> dict:
    return {
        "assay": "eclip",
        "preprocessing": {"input_stage": "raw"},
        "modules": {
            "inline_barcode_demultiplexing": {
                "enabled": True,
                "read": "R1",
                "end": "5prime",
                "max_mismatches": 0,
            },
            "umi_extraction": {
                "enabled": True,
                "method": "string",
                "read1": {"role": "biological"},
                "read2": {
                    "pattern": "NNNNN",
                    "location": "5prime",
                    "role": "biological",
                },
            },
            "umi_deduplication": {
                "enabled": True,
                "method": "unique",
                "edit_distance": 0,
                "umi_tag": "RX",
                "qname_regex": r"_(?P<umi>[ACGTN]+)$",
            },
            "enrichment": {"enabled": False},
        },
    }


def eclip_errors(
    config: dict, rows: list[dict[str, str]]
) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    _validate_eclip(config, rows, errors, warnings)
    return errors, warnings


def test_barcode_list_preserves_multiple_assignments() -> None:
    assert parse_inline_barcodes("A01=AAGCAAT,B06=GGCTTGT") == [
        ("A01", "AAGCAAT"),
        ("B06", "GGCTTGT"),
    ]
    with pytest.raises(ValueError, match="Name=SEQUENCE"):
        parse_inline_barcodes("AAGCAAT")


def test_umi_only_mate_layouts_and_aligner_read_helpers() -> None:
    r2_umi_cfg = {
        "preprocessing": {"input_stage": "raw"},
        "modules": {
            "umi_extraction": {
                "enabled": True,
                "read2": {"pattern": "NNNNNNNN", "role": "umi_only"},
            }
        },
    }
    rows = [sample("ip")]
    assert processed_unit_layout(r2_umi_cfg, rows, "ip__L001") == "single"
    assert processed_sample_layout(r2_umi_cfg, rows, "ip") == "single"

    r1_umi_cfg = {
        "preprocessing": {"input_stage": "raw"},
        "modules": {
            "umi_extraction": {
                "enabled": True,
                "read1": {"pattern": "NNNNNNNN", "role": "umi_only"},
                "read2": {"role": "biological"},
            }
        },
    }
    assert processed_unit_layout(r1_umi_cfg, rows, "ip__L001") == "single_r2"
    assert processed_sample_layout(r1_umi_cfg, rows, "ip") == "single"
    r1_list = alignment_fastq_1(r1_umi_cfg, rows, "ip", "results", True)
    r2_list = alignment_fastq_2(r1_umi_cfg, rows, "ip", "results", True)
    assert r1_list == []
    assert r2_list == ["results/fastq/trimmed/ip__L001_R2.fastq.gz"]
    assert star_reads_arg(r1_list, r2_list) == "results/fastq/trimmed/ip__L001_R2.fastq.gz"
    assert bowtie2_reads_arg(r1_list, r2_list) == "-U results/fastq/trimmed/ip__L001_R2.fastq.gz"
    assert bwa_mem2_reads_arg(r1_list, r2_list) == "<(gzip -cdf results/fastq/trimmed/ip__L001_R2.fastq.gz)"


def test_clip_strand_flags_select_one_mate_and_orientation() -> None:
    assert strand_flag_args("R2", "alignment", "plus", "paired") == (
        "--samFlagInclude 128 --samFlagExclude 16"
    )
    assert strand_flag_args("R2", "alignment", "minus", "paired") == "--samFlagInclude 144"
    assert strand_flag_args("R2", "reverse", "plus", "paired") == "--samFlagInclude 144"
    assert (
        strand_flag_args("R1", "alignment", "plus", "single") == "--samFlagExclude 16"
    )
    assert crosslink_offset({"modules": {"crosslink_sites": {"end": "5prime"}}}) == 1
    assert crosslink_offset({"modules": {"crosslink_sites": {"end": "3prime"}}}) == -1


def test_representative_count_samtools_args_matches_strand_combined_mate() -> None:
    cfg = {
        "preprocessing": {"input_stage": "raw"},
        "tools": {"deeptools": {"params": {"min_mapping_quality": 10}}},
        "modules": {
            "coverage": {
                "read": "R2",
                "params": {"min_mapping_quality": 20, "ignore_duplicates": True},
            }
        },
    }
    args = representative_count_samtools_args(cfg, [sample("ip")], "ip", "coverage")
    assert args == "-F 3844 -f 128 -q 20"


def test_metadata_derived_eclip_adapters_and_dedup_args() -> None:
    cfg = deepcopy(eclip_config())
    cfg["modules"]["trimming"] = {
        "enabled": True,
        "tool": "cutadapt",
        "adapters": {
            "mode": "sequences",
            "barcode_aware": True,
            "read1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "read2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        },
        "params": {"quality_cutoff": 6, "minimum_length": 18, "minimum_overlap": 5, "times": 2},
    }
    rows = [sample("ip", inline_barcodes="A01=AAGCAAT,B06=GGCTTGT")]
    overrides = eclip_barcode_adapter_overrides(cfg, rows, "ip__L001", True)
    assert overrides["adapter_r1"] == "NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=6"
    assert overrides["front_r1"] == ["CTTCCGATCTAAGCAAT", "CTTCCGATCTGGCTTGT"]
    # revcomp(AAGCAAT)=ATTGCTT, revcomp(GGCTTGT)=ACAAGCC
    assert overrides["adapter_r2"][0] == "ATTGCTTAGATCGGA"
    assert any(chunk.startswith("ACAAGCCAGATCGGA") for chunk in overrides["adapter_r2"])
    rendered = render_tool_params(
        cfg,
        "cutadapt",
        overrides=trimming_adapter_params(cfg, "cutadapt", "paired", rows, "ip__L001"),
        module="trimming",
    )
    assert "--match-read-wildcards" in rendered
    assert "-O 5" in rendered
    assert "-q 6" in rendered

    unique_dedup = dedup_args(cfg, paired=True, per_barcode=True)
    assert "--method=unique" in unique_dedup
    assert "--edit-distance-threshold" not in unique_dedup
    assert "--per-cell --cell-tag=BC" in unique_dedup

    cfg["modules"]["umi_deduplication"]["method"] = "directional"
    cfg["modules"]["umi_deduplication"]["edit_distance"] = 1
    dir_dedup = dedup_args(cfg, paired=False, per_barcode=False)
    assert "--method=directional" in dir_dedup
    assert "--edit-distance-threshold=1" in dir_dedup
    assert "--per-cell" not in dir_dedup


def test_stamp_fastq_and_qname_umi_to_rx_support_both_qname_anchors(
    tmp_path: Path,
) -> None:
    fq_in = tmp_path / "in.fastq.gz"
    fq_out = tmp_path / "stamped.fastq.gz"
    with gzip.open(fq_in, "wt", encoding="ascii") as handle:
        handle.write("@AAAAA:read_one 1:N:0:1\nACGTACGT\n+\nIIIIIIII\n")

    subprocess.run(
        [
            sys.executable,
            str(ROOT / "workflow" / "scripts" / "stamp_fastq_barcodes.py"),
            "--entry",
            f"A01=AAGCAAT={fq_in}",
            "--output",
            str(fq_out),
        ],
        check=True,
    )
    with gzip.open(fq_out, "rt", encoding="ascii") as handle:
        header = handle.readline().strip()
    assert header == "@AAAAA:read_one|BC=A01:AAGCAAT 1:N:0:1"

    bam_in = tmp_path / "in.bam"
    bam_out = tmp_path / "out.bam"
    sam_header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}
    with pysam.AlignmentFile(bam_in, "wb", header=sam_header) as out_bam:
        rec = pysam.AlignedSegment()
        rec.query_name = "AAAAA:read_one|BC=A01:AAGCAAT"
        rec.query_sequence = "ACGTACGT"
        rec.flag = 0
        rec.reference_id = 0
        rec.reference_start = 100
        rec.mapping_quality = 255
        rec.cigar = ((0, 8),)
        rec.query_qualities = pysam.qualitystring_to_array("IIIIIIII")
        out_bam.write(rec)

    subprocess.run(
        [
            sys.executable,
            str(ROOT / "workflow" / "scripts" / "qname_umi_to_rx.py"),
            "--input",
            str(bam_in),
            "--output",
            str(bam_out),
            "--regex",
            r"^(?P<umi>[ACGTN]+):",
            "--tag",
            "UR",
        ],
        check=True,
    )
    with pysam.AlignmentFile(bam_out, "rb") as in_bam:
        tagged = next(in_bam.fetch(until_eof=True))
    assert tagged.query_name == "AAAAA:read_one|BC=A01"
    assert tagged.get_tag("UR") == "AAAAA"
    assert tagged.get_tag("BC") == "AAGCAAT"


def test_poly_a_suppresses_none_adapter_warning_and_rejects_fastp() -> None:
    cfg = {
        "tools": {"cutadapt": {"params": {}}},
        "modules": {
            "trimming": {
                "enabled": True,
                "tool": "cutadapt",
                "adapters": {"mode": "none"},
                "params": {"poly_a": True},
            }
        },
    }
    errors: list[str] = []
    warnings: list[str] = []
    _validate_adapters(cfg, [sample("ip")], errors, warnings)
    assert errors == []
    assert warnings == []
    assert "--poly-a" in render_tool_params(cfg, "cutadapt", module="trimming")

    cfg["modules"]["trimming"]["tool"] = "fastp"
    _validate_adapters(cfg, [sample("ip")], errors, warnings)
    assert any("poly_a is supported only" in err for err in errors)


def test_multistage_scalar_extra_is_rejected() -> None:
    cfg = {"tools": {"cutadapt": {"extra": "--anywhere"}}}
    errors: list[str] = []
    _validate_tool_extras(cfg, errors)
    assert any("must be a stage mapping" in err for err in errors)
    with pytest.raises(ValueError, match="multiple stages"):
        tool_extra(cfg, "cutadapt", "trim")


def test_eclip_sra_requires_explicit_sra_layout_and_preprocessed_barcodes() -> None:
    sra_row = sample("ip")
    sra_row["fastq_1"] = ""
    sra_row["fastq_2"] = ""
    sra_row["sra_id"] = "SRR123456"
    sra_row["sra_layout"] = ""
    errors, _ = eclip_errors(eclip_config(), [sra_row])
    assert any("requires explicit sra_layout" in err for err in errors)

    pre_cfg = deepcopy(eclip_config())
    pre_cfg["preprocessing"]["input_stage"] = "umi_extracted"
    errors_missing, _ = eclip_errors(pre_cfg, [sample("ip", inline_barcodes="")])
    assert any("requires one known Name=SEQUENCE" in err for err in errors_missing)

    errors_multi, _ = eclip_errors(
        pre_cfg, [sample("ip", inline_barcodes="A01=AAGCAAT,B06=GGCTTGT")]
    )
    assert any("at most one Name=SEQUENCE" in err for err in errors_multi)


def test_only_sm_input_is_an_enrichment_denominator() -> None:
    rows = [
        sample("ip", control="input"),
        sample("input", library_type="sm_input"),
        sample("tag", library_type="tag_control"),
    ]
    assert enrichment_pairs(rows) == [("ip", "input")]


def test_module_params_override_base_tool_profile() -> None:
    config = {
        "tools": {"deeptools": {"params": {"bin_size": 25, "normalize_using": "CPM"}}},
        "modules": {"coverage": {"params": {"bin_size": 10}}},
    }
    rendered = render_tool_params(config, "deeptools", module="coverage")
    assert "--binSize 10" in rendered
    assert "--normalizeUsing CPM" in rendered


def test_qname_regex_requires_named_umi_group() -> None:
    validate_qname_regex(r"_(?P<umi>[ACGTN]+)$")
    with pytest.raises(ValueError, match="named group"):
        validate_qname_regex(r"_([ACGTN]+)$")


def test_shared_pool_requires_identical_sources() -> None:
    rows = [sample("ip_a"), sample("ip_b")]
    for row, barcode in zip(rows, ("BC_A=AAAAAAA", "BC_B=CCCCCCC"), strict=True):
        row["pool_id"] = "pool_1"
        row["inline_barcodes"] = barcode
    rows[1]["fastq_1"] = "different.fastq.gz"
    errors, _ = eclip_errors(eclip_config(), rows)
    assert any("repeat identical FASTQ/SRA inputs" in error for error in errors)


def test_duplicate_barcode_sequence_in_pool_is_rejected() -> None:
    rows = [sample("ip_a"), sample("ip_b")]
    for index, row in enumerate(rows):
        row["pool_id"] = "pool_1"
        row["inline_barcodes"] = f"BC{index}=AAGCAAT"
    errors, _ = eclip_errors(eclip_config(), rows)
    assert any("sequence AAGCAAT has multiple assignments" in error for error in errors)


def test_tag_control_cannot_be_an_enrichment_denominator() -> None:
    rows = [sample("ip", control="tag"), sample("tag", library_type="tag_control")]
    errors, _ = eclip_errors(eclip_config(), rows)
    assert any("must reference an sm_input" in error for error in errors)


def test_regex_extraction_requires_umi_tools_named_group() -> None:
    config = deepcopy(eclip_config())
    config["modules"]["umi_extraction"] = {
        "enabled": True,
        "method": "regex",
        "read1": {
            "pattern": r"(?P<discard_1>AAA)(?P<bases>[ACGT]+)",
            "role": "biological",
        },
    }
    errors, _ = eclip_errors(config, [sample("ip")])
    assert any("requires a named UMI group" in error for error in errors)


def test_unsupported_unstranded_umi_deduplication_is_rejected() -> None:
    config = deepcopy(eclip_config())
    config["modules"]["umi_deduplication"]["strand_aware"] = False
    errors, _ = eclip_errors(config, [sample("ip")])
    assert any("strand_aware must be true" in error for error in errors)
