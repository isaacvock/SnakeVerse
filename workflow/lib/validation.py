from __future__ import annotations

from pathlib import Path
from typing import Any

from refs import INDEX_KEYS, genome_fasta, index_requires_build, path_exists
from samples import (
    featurecounts_paired_end,
    featurecounts_strand,
    has_sra_samples,
    required_columns_for_assay,
    sample_ids,
    sample_layout,
)


TRIMMING_TOOLS = {"cutadapt", "fastp"}
ALIGNMENT_TOOLS = set(INDEX_KEYS)
FEATURECOUNTS_RUN_KEYS = {
    "feature_type",
    "attribute_type",
    "strand",
    "strandedness",
    "paired_end",
    "count_read_pairs",
    "require_both_ends_mapped",
    "count_multimapping_reads",
    "count_overlapping_features",
}
ADAPTER_TOOL_KEYS = {
    "fastp": {"detect_adapter_for_pe", "adapter_sequence", "adapter_sequence_r2"},
    "cutadapt": {"adapter_r1", "adapter_r2"},
}


def _module(config: dict[str, Any], name: str) -> dict[str, Any]:
    value = (config.get("modules") or {}).get(name, {})
    return value if isinstance(value, dict) else {}


def _enabled(config: dict[str, Any], name: str) -> bool:
    return bool(_module(config, name).get("enabled", False))


def _project_path(config: dict[str, Any], value: str) -> Path:
    root = Path(config.get("_ngsflow", {}).get("project_root", "."))
    path = Path(value)
    return path if path.is_absolute() else root / path


def _required_tools(config: dict[str, Any], samples: list[dict[str, str]]) -> list[str]:
    tools: list[str] = []
    if _enabled(config, "fastq_qc"):
        tools.append("fastqc")
    if _enabled(config, "sra_download") and has_sra_samples(samples):
        tools.append("sra_tools")
    if _enabled(config, "trimming"):
        tools.append(str(_module(config, "trimming").get("tool") or ""))
    if _enabled(config, "alignment"):
        tools.extend([str(_module(config, "alignment").get("tool") or ""), "samtools"])
    if any(_enabled(config, name) for name in ("mark_duplicates", "bam_filter", "bam_qc")):
        tools.append("samtools")
    if _enabled(config, "coverage"):
        tools.append("deeptools")
    if _enabled(config, "multiqc"):
        tools.append("multiqc")
    if _enabled(config, "gene_counts"):
        tools.append(str(_module(config, "gene_counts").get("tool") or "featurecounts"))
    if _enabled(config, "salmon_quantification"):
        tools.append("salmon")
    if _enabled(config, "rsem_quantification"):
        tools.append("rsem")
    if _enabled(config, "peak_calling"):
        tools.append(str(_module(config, "peak_calling").get("tool") or "macs3"))
    if _enabled(config, "frip") or _enabled(config, "tss_enrichment"):
        tools.append("bedtools")
    return sorted({tool for tool in tools if tool})


def _validate_adapters(
    config: dict[str, Any], samples: list[dict[str, str]], errors: list[str], warnings: list[str]
) -> None:
    if not _enabled(config, "trimming"):
        return
    trimming = _module(config, "trimming")
    tool = str(trimming.get("tool") or "")
    if tool not in TRIMMING_TOOLS:
        errors.append(f"Unsupported modules.trimming.tool: {tool}")
        return
    adapters = trimming.get("adapters", {}) or {}
    mode = adapters.get("mode")
    if mode not in {"auto", "sequences", "none"}:
        errors.append("modules.trimming.adapters.mode must be auto, sequences, or none")
    elif mode == "auto" and tool != "fastp":
        errors.append("Adapter auto-detection is supported only when modules.trimming.tool is fastp")
    elif mode == "sequences":
        if not adapters.get("read1"):
            errors.append("Adapter mode sequences requires modules.trimming.adapters.read1")
        paired = any(sample_layout(samples, sample) == "paired" for sample in sample_ids(samples))
        if paired and not adapters.get("read2"):
            errors.append("Paired-end adapter mode sequences requires modules.trimming.adapters.read2")
    elif mode == "none":
        warnings.append("Trimming is enabled without adapter removal (adapters.mode: none)")

    params = config.get("tools", {}).get(tool, {}).get("params", {}) or {}
    conflicts = sorted(ADAPTER_TOOL_KEYS.get(tool, set()).intersection(params))
    if conflicts:
        errors.append(
            f"Adapter settings belong in modules.trimming.adapters; remove from tools/{tool}.yaml: "
            + ", ".join(conflicts)
        )


def validate_resolved_config(
    config: dict[str, Any], samples: list[dict[str, str]]
) -> tuple[list[str], list[str]]:
    """Return contextual validation errors and warnings for a v2 run."""
    errors: list[str] = []
    warnings: list[str] = []
    assay = config.get("assay")

    if config.get("schema_version") != 2:
        errors.append("Run config must set schema_version: 2")
    if assay not in {"generic", "rnaseq", "atacseq"}:
        errors.append(f"Unsupported assay: {assay}")
    for key in ("project", "samples", "reference", "results_dir", "modules"):
        if key not in config:
            errors.append(f"Resolved config is missing required key: {key}")

    if not _enabled(config, "alignment"):
        errors.append("modules.alignment must be enabled")
    aligner = str(_module(config, "alignment").get("tool") or "")
    if aligner not in ALIGNMENT_TOOLS:
        errors.append(f"Unsupported modules.alignment.tool: {aligner}")

    _validate_adapters(config, samples, errors, warnings)

    active_tools = config.get("tools", {}) or {}
    missing_tools = [tool for tool in _required_tools(config, samples) if tool not in active_tools]
    if missing_tools:
        errors.append(
            "Missing active tool configurations: "
            + ", ".join(missing_tools)
            + ". Add each with: python config/bin/ngsflow.py add-tool <tool>"
        )

    reference = config.get("reference", {}) or {}
    if aligner in INDEX_KEYS and index_requires_build(config, aligner) and not genome_fasta(config):
        errors.append(
            f"reference.fasta is required because the {aligner} index must be built"
        )

    if assay == "rnaseq" and _enabled(config, "gene_counts"):
        modes = _module(config, "gene_counts").get("modes", {}) or {}
        if not any(bool(value) for value in modes.values()):
            errors.append("modules.gene_counts is enabled but no counting mode is enabled")
        if not reference.get("gtf"):
            errors.append("reference.gtf is required by modules.gene_counts")
        if _module(config, "gene_counts").get("tool", "featurecounts") != "featurecounts":
            errors.append("modules.gene_counts currently supports only featurecounts")
        try:
            featurecounts_paired_end(samples, config)
            featurecounts_strand(samples, config)
        except ValueError as exc:
            errors.append(str(exc))
        fc_params = active_tools.get("featurecounts", {}).get("params", {}) or {}
        conflicts = sorted(FEATURECOUNTS_RUN_KEYS.intersection(fc_params))
        if conflicts:
            errors.append(
                "Counting semantics belong in modules.gene_counts; remove from "
                "tools/featurecounts.yaml: " + ", ".join(conflicts)
            )

    transcript_quant = _enabled(config, "salmon_quantification") or _enabled(
        config, "rsem_quantification"
    )
    for name in ("salmon_quantification", "rsem_quantification"):
        settings = _module(config, name)
        if _enabled(config, name) and not (
            settings.get("gene_results", True) or settings.get("isoform_results", True)
        ):
            errors.append(f"modules.{name} is enabled but requests no result type")
    if assay == "rnaseq" and transcript_quant:
        if not reference.get("gtf"):
            errors.append("reference.gtf is required by transcript quantification")
        if not reference.get("fasta"):
            errors.append("reference.fasta is required by transcript quantification")

    needs_transcriptome = bool(_module(config, "alignment").get("transcriptome_bam")) or transcript_quant
    if needs_transcriptome:
        if aligner != "star":
            errors.append("STAR is required for transcriptome BAM output")
        star_params = active_tools.get("star", {}).get("params", {}) or {}
        align_params = star_params.get("align", star_params)
        if "TranscriptomeSAM" not in str(align_params.get("quantMode", "")):
            errors.append(
                "STAR transcriptome BAM output requires quantMode to include TranscriptomeSAM"
            )

    if assay == "atacseq":
        replicate_values = {row.get("replicate") for row in samples if row.get("replicate")}
        if len(replicate_values) < 2:
            warnings.append("ENCODE ATAC-seq standards recommend two or more biological replicates")
        if _enabled(config, "bam_filter") and not reference.get("blacklist"):
            warnings.append(
                "ATAC-seq blacklist filtering is recommended; reference.blacklist is blank"
            )
        if _enabled(config, "tss_enrichment"):
            if not reference.get("tss_bed"):
                errors.append("reference.tss_bed is required by modules.tss_enrichment")
            if not reference.get("chrom_sizes"):
                errors.append("reference.chrom_sizes is required by modules.tss_enrichment")
        if _enabled(config, "frip") and not _enabled(config, "peak_calling"):
            errors.append("modules.frip requires modules.peak_calling")

    for ref_key in ("fasta", "gtf", "chrom_sizes", "blacklist", "tss_bed"):
        value = reference.get(ref_key)
        if value and not path_exists(config.get("_ngsflow", {}).get("project_root", "."), str(value)):
            warnings.append(f"Reference path does not exist yet: reference.{ref_key} = {value}")

    required_cols = required_columns_for_assay(str(assay))
    for sample in sample_ids(samples):
        try:
            sample_layout(samples, sample)
        except ValueError as exc:
            errors.append(str(exc))
    for row in samples:
        sample = row.get("sample_id", "<unknown>")
        for column in required_cols:
            if column == "fastq_1" and row.get("sra_id"):
                continue
            if not row.get(column):
                errors.append(f"Sample {sample} is missing {column}")
        if assay == "rnaseq" and row.get("strandedness") not in {
            "unstranded",
            "forward",
            "reverse",
        }:
            errors.append(
                f"Sample {sample} has invalid strandedness '{row.get('strandedness')}'"
            )
        if row.get("sra_layout") and row.get("sra_layout") not in {"single", "paired"}:
            errors.append(f"Sample {sample} has invalid sra_layout '{row.get('sra_layout')}'")
        if row.get("sra_id") and not _enabled(config, "sra_download"):
            errors.append(f"Sample {sample} uses sra_id but modules.sra_download is disabled")
        for fastq_col in ("fastq_1", "fastq_2"):
            fastq = row.get(fastq_col)
            if fastq and not _project_path(config, fastq).exists():
                warnings.append(f"FASTQ path does not exist yet: {fastq}")

    return errors, warnings
