from __future__ import annotations

import sys
from copy import deepcopy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib"))

from clip import crosslink_offset, strand_flag_args, validate_qname_regex  # noqa: E402
from params import render_tool_params  # noqa: E402
from samples import (  # noqa: E402
    enrichment_pairs,
    parse_inline_barcodes,
    processed_sample_layout,
)
from validation import _validate_eclip  # noqa: E402


def sample(
    sample_id: str, library_type: str = "ip", control: str = ""
) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "unit_id": "L001",
        "fastq_1": "R1.fastq.gz",
        "fastq_2": "R2.fastq.gz",
        "condition": "K562",
        "replicate": "1",
        "library_type": library_type,
        "target": "PUM2",
        "pool_id": "",
        "input_sample_id": control,
        "inline_barcodes": "",
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
    assert parse_inline_barcodes("A01=ATTGCTT,B06=ACAAGCC") == [
        ("A01", "ATTGCTT"),
        ("B06", "ACAAGCC"),
    ]
    with pytest.raises(ValueError, match="Name=SEQUENCE"):
        parse_inline_barcodes("ATTGCTT")


def test_umi_only_r2_changes_processed_layout_to_single() -> None:
    config = {
        "preprocessing": {"input_stage": "raw"},
        "modules": {
            "umi_extraction": {
                "enabled": True,
                "read2": {"pattern": "NNNNNNNN", "role": "umi_only"},
            }
        },
    }
    assert processed_sample_layout(config, [sample("ip")], "ip") == "single"


def test_clip_strand_flags_select_one_mate_and_orientation() -> None:
    assert strand_flag_args("R2", "reverse", "plus", "paired") == "--samFlagInclude 144"
    assert strand_flag_args("R2", "reverse", "minus", "paired") == (
        "--samFlagInclude 128 --samFlagExclude 16"
    )
    assert (
        strand_flag_args("R1", "alignment", "plus", "single") == "--samFlagExclude 16"
    )
    assert crosslink_offset({"modules": {"crosslink_sites": {"end": "5prime"}}}) == 1
    assert crosslink_offset({"modules": {"crosslink_sites": {"end": "3prime"}}}) == -1


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
        row["inline_barcodes"] = f"BC{index}=ATTGCTT"
    errors, _ = eclip_errors(eclip_config(), rows)
    assert any("sequence ATTGCTT has multiple assignments" in error for error in errors)


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
