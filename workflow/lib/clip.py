from __future__ import annotations

import re
import shlex
from typing import Any

from samples import (
    parse_inline_barcodes,
    pool_rows,
    processed_sample_layout,
    sample_metadata,
)


def barcode_adapter_args(
    config: dict[str, Any], samples: list[dict[str, str]], pool: str
) -> str:
    module = (config.get("modules", {}) or {}).get(
        "inline_barcode_demultiplexing", {}
    ) or {}
    read = str(module.get("read", "R1")).upper()
    end = str(module.get("end", "5prime"))
    flag = {
        ("R1", "5prime"): "-g",
        ("R1", "3prime"): "-a",
        ("R2", "5prime"): "-G",
        ("R2", "3prime"): "-A",
    }[(read, end)]
    anchor = "^" if end == "5prime" else "$"
    adapters: list[str] = []
    seen: set[str] = set()
    for row in pool_rows(samples, pool):
        for name, sequence in parse_inline_barcodes(row.get("inline_barcodes", "")):
            if name in seen:
                continue
            seen.add(name)
            expression = (
                f"{name}={anchor}{sequence}"
                if end == "5prime"
                else f"{name}={sequence}{anchor}"
            )
            adapters.extend([flag, shlex.quote(expression)])
    mismatches = int(module.get("max_mismatches", 0))
    adapters.extend(["-e", str(mismatches)])
    if not bool(module.get("allow_indels", False)):
        adapters.append("--no-indels")
    return " ".join(adapters)


def demultiplexed_barcode_files(
    results_dir: str, samples: list[dict[str, str]], unit: str, read: str
) -> str:
    from samples import sample_by_unit

    row = sample_by_unit(samples, unit)
    pool = row.get("pool_id", "")
    names = [name for name, _ in parse_inline_barcodes(row.get("inline_barcodes", ""))]
    return " ".join(
        shlex.quote(
            f"{results_dir}/fastq/demultiplexed/pools/{pool}/{name}_{read}.fastq.gz"
        )
        for name in names
    )


def umi_extract_args(config: dict[str, Any], paired: bool) -> str:
    module = (config.get("modules", {}) or {}).get("umi_extraction", {}) or {}
    method = str(module.get("method", "string"))
    parts = ["--extract-method", method]
    locations: set[str] = set()
    for read, flag in (("read1", "--bc-pattern"), ("read2", "--bc-pattern2")):
        settings = module.get(read, {}) or {}
        pattern = settings.get("pattern")
        if not pattern:
            continue
        parts.extend([flag, shlex.quote(str(pattern))])
        locations.add(str(settings.get("location", "5prime")))
    if method == "string" and locations == {"3prime"}:
        parts.append("--3prime")
    return " ".join(parts)


def dedup_args(config: dict[str, Any], paired: bool) -> str:
    module = (config.get("modules", {}) or {}).get("umi_deduplication", {}) or {}
    parts = [
        "--extract-umi-method=tag",
        f"--umi-tag={shlex.quote(str(module.get('umi_tag', 'RX')))}",
        f"--method={shlex.quote(str(module.get('method', 'directional')))}",
        f"--edit-distance-threshold={int(module.get('edit_distance', 1))}",
    ]
    if bool(module.get("spliced_is_unique", True)):
        parts.append("--spliced-is-unique")
    if paired:
        parts.append("--paired")
    return " ".join(parts)


def _mapped_reverse(strand: str, interpretation: str) -> bool:
    reverse = strand == "minus"
    return not reverse if interpretation == "reverse" else reverse


def strand_flag_args(read: str, interpretation: str, strand: str, layout: str) -> str:
    include = 0
    exclude = 0
    if layout == "paired":
        if read == "R1":
            include |= 64
        elif read == "R2":
            include |= 128
    if _mapped_reverse(strand, interpretation):
        include |= 16
    else:
        exclude |= 16
    parts: list[str] = []
    if include:
        parts.extend(["--samFlagInclude", str(include)])
    if exclude:
        parts.extend(["--samFlagExclude", str(exclude)])
    return " ".join(parts)


def track_strand_args(
    config: dict[str, Any],
    samples: list[dict[str, str]],
    sample: str,
    strand: str,
    module: str,
) -> str:
    settings = (config.get("modules", {}) or {}).get(module, {}) or {}
    read = str(settings.get("read", "any"))
    interpretation = str(settings.get("strand", "alignment"))
    return strand_flag_args(
        read, interpretation, strand, processed_sample_layout(config, samples, sample)
    )


def crosslink_offset(config: dict[str, Any]) -> int:
    settings = (config.get("modules", {}) or {}).get("crosslink_sites", {}) or {}
    return 1 if settings.get("end", "5prime") == "5prime" else -1


def qname_umi_regex(config: dict[str, Any]) -> str:
    settings = (config.get("modules", {}) or {}).get("umi_deduplication", {}) or {}
    return str(settings.get("qname_regex") or r"_(?P<umi>[ACGTN]+)$")


def enrichment_control(samples: list[dict[str, str]], sample: str) -> str:
    return sample_metadata(samples, sample, "input_sample_id")


def validate_qname_regex(value: str) -> None:
    compiled = re.compile(value)
    if "umi" not in compiled.groupindex:
        raise ValueError("UMI QNAME regex must contain a named group '(?P<umi>...)'")
