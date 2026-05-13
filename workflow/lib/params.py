from __future__ import annotations

import shlex
from collections.abc import Mapping
from typing import Any

from samples import featurecounts_count_read_pairs, featurecounts_paired_end


def _dash(name: str) -> str:
    return name.replace("_", "-")


def _camel_deeptools(name: str) -> str:
    special = {
        "bin_size": "binSize",
        "normalize_using": "normalizeUsing",
        "effective_genome_size": "effectiveGenomeSize",
        "min_mapping_quality": "minMappingQuality",
        "ignore_duplicates": "ignoreDuplicates",
        "extend_reads": "extendReads",
        "filter_rna_strand": "filterRNAstrand",
    }
    return special.get(name, _dash(name))


def _append_flag(parts: list[str], flag: str, value: Any) -> None:
    if value is None or value == "":
        return
    if isinstance(value, bool):
        if value:
            parts.append(flag)
        return
    if isinstance(value, (list, tuple)):
        parts.append(flag)
        parts.extend(shlex.quote(str(item)) for item in value)
        return
    parts.extend([flag, shlex.quote(str(value))])


def render_generic_params(params: Mapping[str, Any]) -> str:
    parts: list[str] = []
    for key, value in params.items():
        _append_flag(parts, f"--{_dash(key)}", value)
    return " ".join(parts)


def render_fastqc(params: Mapping[str, Any]) -> str:
    return render_generic_params(params)


def render_cutadapt(params: Mapping[str, Any]) -> str:
    flag_map = {
        "adapter_r1": "-a",
        "adapter_r2": "-A",
        "front_r1": "-g",
        "front_r2": "-G",
        "quality_cutoff": "-q",
        "minimum_length": "-m",
        "error_rate": "-e",
        "times": "-n",
        "nextseq_trim": "--nextseq-trim",
    }
    parts: list[str] = []
    for key, value in params.items():
        flag = flag_map.get(key, f"--{_dash(key)}")
        _append_flag(parts, flag, value)
    return " ".join(parts)


def render_fastp(params: Mapping[str, Any]) -> str:
    flag_map = {
        "qualified_quality_phred": "--qualified_quality_phred",
        "unqualified_percent_limit": "--unqualified_percent_limit",
        "length_required": "--length_required",
        "cut_front": "--cut_front",
        "cut_tail": "--cut_tail",
        "cut_right": "--cut_right",
        "detect_adapter_for_pe": "--detect_adapter_for_pe",
        "adapter_sequence": "--adapter_sequence",
        "adapter_sequence_r2": "--adapter_sequence_r2",
    }
    parts: list[str] = []
    for key, value in params.items():
        flag = flag_map.get(key, f"--{_dash(key)}")
        _append_flag(parts, flag, value)
    return " ".join(parts)


def render_bowtie2(params: Mapping[str, Any]) -> str:
    parts: list[str] = []
    for key, value in params.items():
        if key == "sensitivity":
            if value:
                parts.append(str(value))
        elif key == "max_insert_size":
            _append_flag(parts, "-X", value)
        elif key == "no_mixed":
            _append_flag(parts, "--no-mixed", value)
        elif key == "no_discordant":
            _append_flag(parts, "--no-discordant", value)
        else:
            _append_flag(parts, f"--{_dash(key)}", value)
    return " ".join(parts)


def render_bwa_mem2(params: Mapping[str, Any]) -> str:
    parts: list[str] = []
    flag_map = {
        "min_seed_length": "-k",
        "band_width": "-w",
        "score_match": "-A",
        "mismatch_penalty": "-B",
        "gap_open_penalty": "-O",
        "gap_extension_penalty": "-E",
        "clipping_penalty": "-L",
        "mark_shorter_splits": "-M",
    }
    for key, value in params.items():
        flag = flag_map.get(key, f"--{_dash(key)}")
        _append_flag(parts, flag, value)
    return " ".join(parts)


def render_star(params: Mapping[str, Any]) -> str:
    parts: list[str] = []
    for key, value in params.items():
        if value is None or value == "":
            continue
        flag = f"--{key}"
        if isinstance(value, bool):
            if value:
                parts.append(flag)
            continue
        if isinstance(value, (list, tuple)):
            parts.append(flag)
            parts.extend(shlex.quote(str(item)) for item in value)
            continue
        if isinstance(value, str) and " " in value:
            parts.append(flag)
            parts.extend(shlex.quote(piece) for piece in value.split())
            continue
        parts.extend([flag, shlex.quote(str(value))])
    return " ".join(parts)


def render_samtools_view(params: Mapping[str, Any]) -> str:
    parts: list[str] = []
    if "min_mapq" in params:
        _append_flag(parts, "-q", params["min_mapq"])
    if "required_flags" in params:
        _append_flag(parts, "-f", params["required_flags"])
    if "excluded_flags" in params:
        _append_flag(parts, "-F", params["excluded_flags"])
    if params.get("keep_duplicates") is False and "excluded_flags" not in params:
        _append_flag(parts, "-F", 1024)
    return " ".join(parts)


def render_samtools_sort(params: Mapping[str, Any]) -> str:
    parts: list[str] = []
    if "memory_per_thread" in params:
        _append_flag(parts, "-m", params["memory_per_thread"])
    return " ".join(parts)


def render_featurecounts(params: Mapping[str, Any]) -> str:
    parts: list[str] = []
    flag_map = {
        "annotation_format": "-F",
        "feature_type": "-t",
        "attribute_type": "-g",
        "strand": "-s",
        "count_multimapping_reads": "-M",
        "count_overlapping_features": "-O",
        "paired_end": "-p",
        "count_read_pairs": "--countReadPairs",
        "require_both_ends_mapped": "-B",
        "count_chimeric_fragments": "-C",
        "minimum_overlap": "--minOverlap",
        "non_overlap": "--nonOverlap",
    }
    for key, value in params.items():
        flag = flag_map.get(key, f"--{_dash(key)}")
        _append_flag(parts, flag, value)
    return " ".join(parts)


def render_featurecounts_for_config(
    config: Mapping[str, Any],
    samples: list[dict[str, str]],
    overrides: Mapping[str, Any] | None = None,
    drop_keys: tuple[str, ...] = (),
) -> str:
    params = dict(config.get("tools", {}).get("featurecounts", {}).get("params", {}) or {})
    for key in drop_keys:
        params.pop(key, None)
    params.update(
        {
            "paired_end": featurecounts_paired_end(samples, config),
            "count_read_pairs": featurecounts_count_read_pairs(samples, config),
        }
    )
    if overrides:
        params.update(overrides)
    return render_featurecounts(params)


def render_salmon(params: Mapping[str, Any]) -> str:
    flag_map = {
        "library_type": "-l",
        "seq_bias": "--seqBias",
        "gc_bias": "--gcBias",
        "validate_mappings": "--validateMappings",
        "num_bootstraps": "--numBootstraps",
        "num_gibbs_samples": "--numGibbsSamples",
        "gene_map": "--geneMap",
    }
    parts: list[str] = []
    for key, value in params.items():
        flag = flag_map.get(key, f"--{_dash(key)}")
        _append_flag(parts, flag, value)
    return " ".join(parts)


def render_rsem(params: Mapping[str, Any]) -> str:
    return render_generic_params(params)


def render_deeptools(params: Mapping[str, Any]) -> str:
    parts: list[str] = []
    for key, value in params.items():
        _append_flag(parts, f"--{_camel_deeptools(key)}", value)
    return " ".join(parts)


def render_multiqc(params: Mapping[str, Any]) -> str:
    return render_generic_params(params)


def render_macs3(params: Mapping[str, Any]) -> str:
    parts: list[str] = []
    flag_map = {
        "genome_size": "-g",
        "qvalue": "-q",
        "pvalue": "-p",
        "keep_dup": "--keep-dup",
        "nomodel": "--nomodel",
        "shift": "--shift",
        "extsize": "--extsize",
        "bdg": "-B",
        "trackline": "--trackline",
        "broad": "--broad",
        "broad_cutoff": "--broad-cutoff",
        "cutoff_analysis": "--cutoff-analysis",
        "scale_to": "--scale-to",
        "call_summits": "--call-summits",
        "nolambda": "--nolambda",
    }
    for key, value in params.items():
        flag = flag_map.get(key, f"--{_dash(key)}")
        _append_flag(parts, flag, value)
    return " ".join(parts)


RENDERERS = {
    "fastqc": render_fastqc,
    "cutadapt": render_cutadapt,
    "fastp": render_fastp,
    "bowtie2": render_bowtie2,
    "bwa_mem2": render_bwa_mem2,
    "star": render_star,
    "featurecounts": render_featurecounts,
    "salmon": render_salmon,
    "rsem": render_rsem,
    "deeptools": render_deeptools,
    "macs3": render_macs3,
    "multiqc": render_multiqc,
}


def tool_profile(config: Mapping[str, Any], tool: str) -> Mapping[str, Any]:
    return config.get("tools", {}).get(tool, {})


def tool_extra(config: Mapping[str, Any], tool: str) -> str:
    return str(tool_profile(config, tool).get("extra") or "")


def render_tool_params(
    config: Mapping[str, Any],
    tool: str,
    section: str | None = None,
    overrides: Mapping[str, Any] | None = None,
) -> str:
    profile = tool_profile(config, tool)
    params: Mapping[str, Any] = profile.get("params", {}) or {}
    if section:
        if isinstance(params.get(section), Mapping):
            params = params.get(section, {}) or {}
        elif section == "align":
            params = {
                key: value
                for key, value in params.items()
                if key not in {"align", "index", "sort", "filter"}
            }
        else:
            params = {}
    merged = dict(params)
    if overrides:
        merged.update({key: value for key, value in overrides.items() if value is not None})
    if tool == "samtools" and section == "filter":
        return render_samtools_view(merged)
    if tool == "samtools" and section == "sort":
        return render_samtools_sort(merged)
    renderer = RENDERERS.get(tool, render_generic_params)
    return renderer(merged)


def resource_value(config: Mapping[str, Any], rule_name: str, key: str, default: Any) -> Any:
    return config.get("resources", {}).get(rule_name, {}).get(key, default)
