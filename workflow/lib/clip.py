from __future__ import annotations

import re
import shlex
from collections.abc import Mapping
from typing import Any

from samples import (
    parse_inline_barcodes,
    pool_rows,
    processed_sample_layout,
    sample_by_unit,
    sample_metadata,
)

IUPAC_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def reverse_complement_iupac(sequence: str) -> str:
    return sequence.upper().translate(IUPAC_COMPLEMENT)[::-1]


def barcode_adapter_args(
    config: Mapping[str, Any], samples: list[dict[str, str]], pool: str
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


def demultiplexed_barcode_entries(
    results_dir: str, samples: list[dict[str, str]], unit: str, read: str
) -> str:
    row = sample_by_unit(samples, unit)
    pool = row.get("pool_id", "")
    barcodes = parse_inline_barcodes(row.get("inline_barcodes", ""))
    return " ".join(
        "--entry "
        + shlex.quote(
            f"{name}={seq}={results_dir}/fastq/demultiplexed/pools/{pool}/{name}_{read}.fastq.gz"
        )
        for name, seq in barcodes
    )


def single_unit_barcode_entry(
    samples: list[dict[str, str]], unit: str, fastq_path: str
) -> str:
    row = sample_by_unit(samples, unit)
    barcodes = parse_inline_barcodes(row.get("inline_barcodes", ""))
    if len(barcodes) != 1:
        raise ValueError(
            f"Preprocessed unit {unit} must specify exactly one inline barcode, got {len(barcodes)}"
        )
    name, seq = barcodes[0]
    return f"--entry {shlex.quote(f'{name}={seq}={fastq_path}')}"


def eclip_barcode_adapter_overrides(
    config: Mapping[str, Any],
    samples: list[dict[str, str]],
    unit: str,
    paired: bool,
) -> dict[str, Any]:
    """Derive ENCODE eCLIP R1/R2 adapters from unit inline-barcode and UMI metadata."""
    trimming = (config.get("modules", {}) or {}).get("trimming", {}) or {}
    adapters = trimming.get("adapters", {}) or {}
    read1_raw = adapters.get("read1")
    read2_raw = adapters.get("read2")
    read1_base = str(read1_raw) if isinstance(read1_raw, str) else ""
    read2_base = str(read2_raw) if isinstance(read2_raw, str) else ""

    row = sample_by_unit(samples, unit)
    barcodes = parse_inline_barcodes(row.get("inline_barcodes", ""))

    umi_module = (config.get("modules", {}) or {}).get("umi_extraction", {}) or {}
    r2_umi = umi_module.get("read2", {}) or {}
    umi_len = 0
    if (
        str(umi_module.get("method", "string")) == "string"
        and str(r2_umi.get("location", "5prime")) == "5prime"
        and str(r2_umi.get("role", "biological")) == "biological"
        and r2_umi.get("pattern")
    ):
        umi_len = len(str(r2_umi["pattern"]))

    overrides: dict[str, Any] = {"match_read_wildcards": True}
    if read1_base:
        if umi_len > 0 and ";" not in read1_base:
            overrides["adapter_r1"] = (
                f"{'N' * umi_len}{read1_base};min_overlap={umi_len + 1}"
            )
        else:
            overrides["adapter_r1"] = read1_base
    elif read1_raw:
        overrides["adapter_r1"] = read1_raw

    demux = (config.get("modules", {}) or {}).get(
        "inline_barcode_demultiplexing", {}
    ) or {}
    if (
        barcodes
        and str(demux.get("read", "R1")).upper() == "R1"
        and str(demux.get("end", "5prime")) == "5prime"
    ):
        overrides["front_r1"] = [f"CTTCCGATCT{seq}" for _, seq in barcodes]

    if paired:
        if barcodes and read2_base:
            adapter_head = read2_base[:20]
            chunks: list[str] = []
            for _, seq in barcodes:
                full = reverse_complement_iupac(seq) + adapter_head
                window = min(15, len(full))
                chunks.extend(
                    full[idx : idx + window]
                    for idx in range(max(1, len(full) - window + 1))
                )
            overrides["adapter_r2"] = list(dict.fromkeys(chunks))
        elif read2_raw:
            overrides["adapter_r2"] = read2_raw

    return overrides


def umi_extract_args(config: Mapping[str, Any], paired: bool) -> str:
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


def dedup_args(
    config: Mapping[str, Any], paired: bool, per_barcode: bool = False
) -> str:
    module = (config.get("modules", {}) or {}).get("umi_deduplication", {}) or {}
    method = str(module.get("method", "directional"))
    parts = [
        "--extract-umi-method=tag",
        f"--umi-tag={shlex.quote(str(module.get('umi_tag', 'RX')))}",
        f"--method={shlex.quote(method)}",
    ]
    if method != "unique":
        parts.append(
            f"--edit-distance-threshold={int(module.get('edit_distance', 1))}"
        )
    if per_barcode:
        parts.extend(["--per-cell", "--cell-tag=BC"])
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
    config: Mapping[str, Any],
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


def representative_count_samtools_args(
    config: Mapping[str, Any],
    samples: list[dict[str, str]],
    sample: str,
    module: str,
) -> str:
    """Return strand-combined samtools view flags matching deepTools' kept reads."""
    settings = (config.get("modules", {}) or {}).get(module, {}) or {}
    deeptools_params = (
        (config.get("tools", {}) or {}).get("deeptools", {}).get("params", {}) or {}
    )
    module_params = settings.get("params", {}) or {}
    if module == "coverage":
        min_mapq = int(
            module_params.get(
                "min_mapping_quality", deeptools_params.get("min_mapping_quality", 0)
            )
            or 0
        )
        ignore_dups = bool(
            module_params.get(
                "ignore_duplicates", deeptools_params.get("ignore_duplicates", False)
            )
        )
    else:
        min_mapq = int(settings.get("min_mapping_quality", 10) or 0)
        ignore_dups = False

    excluded = 2820 | (1024 if ignore_dups else 0)
    parts = ["-F", str(excluded)]
    layout = processed_sample_layout(config, samples, sample)
    read = str(settings.get("read", "any"))
    if layout == "paired":
        if read == "R1":
            parts.extend(["-f", "64"])
        elif read == "R2":
            parts.extend(["-f", "128"])
    if min_mapq > 0:
        parts.extend(["-q", str(min_mapq)])
    return " ".join(parts)


def coverage_normalization_mode(config: Mapping[str, Any]) -> str:
    deeptools_params = (
        (config.get("tools", {}) or {}).get("deeptools", {}).get("params", {}) or {}
    )
    module_params = (
        ((config.get("modules", {}) or {}).get("coverage", {}) or {}).get("params", {})
        or {}
    )
    value = module_params.get(
        "normalize_using", deeptools_params.get("normalize_using", "CPM")
    )
    return str(value or "None")


def crosslink_offset(config: Mapping[str, Any]) -> int:
    settings = (config.get("modules", {}) or {}).get("crosslink_sites", {}) or {}
    return 1 if settings.get("end", "5prime") == "5prime" else -1


def qname_umi_regex(config: Mapping[str, Any]) -> str:
    settings = (config.get("modules", {}) or {}).get("umi_deduplication", {}) or {}
    return str(settings.get("qname_regex") or r"_(?P<umi>[ACGTN]+)$")


def umi_sam_tag(config: Mapping[str, Any]) -> str:
    settings = (config.get("modules", {}) or {}).get("umi_deduplication", {}) or {}
    return str(settings.get("umi_tag") or "RX")


def enrichment_control(samples: list[dict[str, str]], sample: str) -> str:
    return sample_metadata(samples, sample, "input_sample_id")


def validate_qname_regex(value: str) -> None:
    compiled = re.compile(value)
    if "umi" not in compiled.groupindex:
        raise ValueError("UMI QNAME regex must contain a named group '(?P<umi>...)'")
