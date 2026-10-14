from __future__ import annotations

from pathlib import Path
import re
from typing import Any

from refs import INDEX_KEYS, genome_fasta, index_requires_build, path_exists
from samples import (
    featurecounts_paired_end,
    featurecounts_strand,
    has_sra_samples,
    required_columns_for_assay,
    input_stage,
    parse_inline_barcodes,
    processed_sample_layout,
    sample_ids,
    sample_metadata,
    sample_layout,
    unit_key,
)
from clip import validate_qname_regex

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
RSEM_RUN_KEYS = {"alignments", "paired_end", "strandedness", "num_threads"}
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
    if has_sra_samples(samples):
        tools.append("sra_tools")
    if _enabled(config, "trimming"):
        tools.append(str(_module(config, "trimming").get("tool") or ""))
    if (
        _enabled(config, "inline_barcode_demultiplexing")
        and input_stage(config) == "raw"
    ):
        tools.append("cutadapt")
    if _enabled(config, "umi_extraction") and input_stage(config) != "umi_extracted":
        tools.append("umi_tools")
    if _enabled(config, "umi_deduplication"):
        tools.extend(["umi_tools", "samtools"])
    if _enabled(config, "alignment"):
        tools.extend([str(_module(config, "alignment").get("tool") or ""), "samtools"])
    if any(
        _enabled(config, name) for name in ("mark_duplicates", "bam_filter", "bam_qc")
    ):
        tools.append("samtools")
    if any(
        _enabled(config, name) for name in ("coverage", "crosslink_sites", "enrichment")
    ):
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
    config: dict[str, Any],
    samples: list[dict[str, str]],
    errors: list[str],
    warnings: list[str],
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
        errors.append(
            "Adapter auto-detection is supported only when modules.trimming.tool is fastp"
        )
    elif mode == "sequences":
        if not adapters.get("read1"):
            errors.append(
                "Adapter mode sequences requires modules.trimming.adapters.read1"
            )
        paired = any(
            processed_sample_layout(config, samples, sample) == "paired"
            for sample in sample_ids(samples)
        )
        if paired and not adapters.get("read2"):
            errors.append(
                "Paired-end adapter mode sequences requires modules.trimming.adapters.read2"
            )
    elif mode == "none":
        warnings.append(
            "Trimming is enabled without adapter removal (adapters.mode: none)"
        )

    params = config.get("tools", {}).get(tool, {}).get("params", {}) or {}
    conflicts = sorted(ADAPTER_TOOL_KEYS.get(tool, set()).intersection(params))
    if conflicts:
        errors.append(
            f"Adapter settings belong in modules.trimming.adapters; remove from tools/{tool}.yaml: "
            + ", ".join(conflicts)
        )


def _validate_eclip(
    config: dict[str, Any],
    samples: list[dict[str, str]],
    errors: list[str],
    warnings: list[str],
) -> None:
    if config.get("assay") != "eclip":
        return

    stage = input_stage(config)
    if stage not in {"raw", "demultiplexed", "umi_extracted"}:
        errors.append(
            "preprocessing.input_stage must be raw, demultiplexed, or umi_extracted"
        )
    if stage != "raw" and _enabled(config, "inline_barcode_demultiplexing"):
        warnings.append(
            f"Inline-barcode demultiplexing is bypassed for preprocessing.input_stage: {stage}"
        )
    if stage == "umi_extracted" and _enabled(config, "umi_extraction"):
        warnings.append(
            "UMI extraction is bypassed for preprocessing.input_stage: umi_extracted"
        )

    if _enabled(config, "mark_duplicates") and _enabled(config, "umi_deduplication"):
        errors.append(
            "modules.mark_duplicates and modules.umi_deduplication are mutually exclusive"
        )
    if (
        _enabled(config, "umi_deduplication")
        and stage != "umi_extracted"
        and not _enabled(config, "umi_extraction")
    ):
        errors.append(
            "UMI deduplication requires UMI extraction unless preprocessing.input_stage is umi_extracted"
        )

    demux = _module(config, "inline_barcode_demultiplexing")
    if _enabled(config, "inline_barcode_demultiplexing"):
        if str(demux.get("read", "R1")) not in {"R1", "R2"}:
            errors.append("modules.inline_barcode_demultiplexing.read must be R1 or R2")
        if str(demux.get("end", "5prime")) not in {"5prime", "3prime"}:
            errors.append(
                "modules.inline_barcode_demultiplexing.end must be 5prime or 3prime"
            )
        try:
            if int(demux.get("max_mismatches", 0)) < 0:
                raise ValueError
        except (TypeError, ValueError):
            errors.append(
                "modules.inline_barcode_demultiplexing.max_mismatches must be >= 0"
            )

    umi = _module(config, "umi_extraction")
    if _enabled(config, "umi_extraction"):
        method = str(umi.get("method", "string"))
        if method not in {"string", "regex"}:
            errors.append("modules.umi_extraction.method must be string or regex")
        configured = 0
        locations: set[str] = set()
        for read in ("read1", "read2"):
            settings = umi.get(read, {}) or {}
            if not isinstance(settings, dict):
                errors.append(f"modules.umi_extraction.{read} must be a mapping")
                continue
            pattern = str(settings.get("pattern") or "")
            if pattern:
                configured += 1
                if method == "string" and not re.fullmatch(r"[NCX]+", pattern):
                    errors.append(
                        f"modules.umi_extraction.{read}.pattern must contain only N, C, or X for string extraction"
                    )
                if method == "regex":
                    try:
                        compiled = re.compile(pattern)
                        if not any(
                            name.startswith("umi_") for name in compiled.groupindex
                        ):
                            errors.append(
                                f"modules.umi_extraction.{read}.pattern requires a named "
                                "UMI group such as (?P<umi_1>...)"
                            )
                    except re.error as exc:
                        errors.append(f"Invalid UMI regex for {read}: {exc}")
                location = str(settings.get("location", "5prime"))
                locations.add(location)
                if location not in {"5prime", "3prime"}:
                    errors.append(
                        f"modules.umi_extraction.{read}.location must be 5prime or 3prime"
                    )
            role = str(settings.get("role", "biological"))
            if role not in {"biological", "umi_only"}:
                errors.append(
                    f"modules.umi_extraction.{read}.role must be biological or umi_only"
                )
        if configured == 0:
            errors.append(
                "modules.umi_extraction requires a read1 and/or read2 pattern"
            )
        if method == "string" and len(locations) > 1:
            errors.append(
                "String UMI extraction requires all configured UMIs at the same read end"
            )
        if str((umi.get("read1", {}) or {}).get("role", "biological")) == "umi_only":
            errors.append("UMI-only R1 is not supported; use R2 as the UMI-only mate")
        if str(
            (umi.get("read2", {}) or {}).get("role", "biological")
        ) == "umi_only" and not ((umi.get("read2", {}) or {}).get("pattern")):
            errors.append(
                "An UMI-only R2 requires modules.umi_extraction.read2.pattern"
            )

    if _enabled(config, "umi_deduplication"):
        dedup = _module(config, "umi_deduplication")
        if str(dedup.get("method", "directional")) not in {
            "unique",
            "percentile",
            "cluster",
            "adjacency",
            "directional",
        }:
            errors.append("Unsupported modules.umi_deduplication.method")
        if dedup.get("strand_aware", True) is not True:
            errors.append(
                "modules.umi_deduplication.strand_aware must be true; "
                "UMI-tools deduplication is strand-aware and does not provide an unstranded mode"
            )
        try:
            if int(dedup.get("edit_distance", 1)) < 0:
                raise ValueError
        except (TypeError, ValueError):
            errors.append("modules.umi_deduplication.edit_distance must be >= 0")
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]", str(dedup.get("umi_tag", "RX"))):
            errors.append(
                "modules.umi_deduplication.umi_tag must be a two-character SAM tag"
            )
        try:
            validate_qname_regex(
                str(dedup.get("qname_regex") or r"_(?P<umi>[ACGTN]+)$")
            )
        except (ValueError, re.error) as exc:
            errors.append(str(exc))

    raw_layouts = {sample_layout(samples, sample) for sample in sample_ids(samples)}
    if len(raw_layouts) > 1:
        errors.append("eCLIP runs must use one raw FASTQ layout")
    for sample in sample_ids(samples):
        try:
            processed_sample_layout(config, samples, sample)
        except ValueError as exc:
            errors.append(str(exc))
        for key in (
            "condition",
            "replicate",
            "library_type",
            "target",
            "input_sample_id",
        ):
            try:
                sample_metadata(samples, sample, key)
            except ValueError as exc:
                errors.append(str(exc))

    known = set(sample_ids(samples))
    for sample in sample_ids(samples):
        try:
            library_type = sample_metadata(samples, sample, "library_type")
            control = sample_metadata(samples, sample, "input_sample_id")
        except ValueError:
            continue
        if library_type not in {"ip", "sm_input", "tag_control"}:
            errors.append(f"Sample {sample} has invalid library_type '{library_type}'")
            continue
        if library_type != "ip" and control:
            errors.append(f"Only IP samples may set input_sample_id: {sample}")
        if library_type == "ip" and not control and _enabled(config, "enrichment"):
            warnings.append(
                f"IP sample {sample} has no SMInput link; enrichment will be skipped"
            )
        if not control:
            continue
        if control == sample:
            errors.append(f"IP sample {sample} cannot use itself as input")
        elif control not in known:
            errors.append(
                f"IP sample {sample} references unknown input sample {control}"
            )
        else:
            try:
                if sample_metadata(samples, control, "library_type") != "sm_input":
                    errors.append(
                        f"IP sample {sample} input_sample_id must reference an sm_input"
                    )
                for key in ("condition", "replicate"):
                    if sample_metadata(samples, sample, key) != sample_metadata(
                        samples, control, key
                    ):
                        warnings.append(
                            f"IP sample {sample} and SMInput {control} have different {key}; honoring explicit link"
                        )
            except ValueError as exc:
                errors.append(str(exc))

    pools: dict[str, list[dict[str, str]]] = {}
    barcode_name_owners: dict[tuple[str, str], str] = {}
    barcode_sequence_owners: dict[tuple[str, str], tuple[str, str]] = {}
    for row in samples:
        pool = row.get("pool_id", "")
        try:
            barcodes = parse_inline_barcodes(row.get("inline_barcodes", ""))
        except ValueError as exc:
            errors.append(f"Sample unit {unit_key(row)}: {exc}")
            barcodes = []
        if (
            stage == "raw"
            and _enabled(config, "inline_barcode_demultiplexing")
            and barcodes
            and not pool
        ):
            errors.append(
                f"Sample unit {unit_key(row)} has inline barcodes but no pool_id"
            )
        if pool:
            pools.setdefault(pool, []).append(row)
        for name, sequence in barcodes:
            owner_key = (pool, name)
            if owner_key in barcode_name_owners and barcode_name_owners[
                owner_key
            ] != unit_key(row):
                errors.append(
                    f"Inline barcode {name} is assigned to multiple units in pool {pool}"
                )
            barcode_name_owners[owner_key] = unit_key(row)
            sequence_key = (pool, sequence)
            previous = barcode_sequence_owners.get(sequence_key)
            assignment = (name, unit_key(row))
            if previous is not None and previous != assignment:
                errors.append(
                    f"Inline barcode sequence {sequence} has multiple assignments in pool {pool}: "
                    f"{previous[0]} ({previous[1]}) and {name} ({unit_key(row)})"
                )
            barcode_sequence_owners[sequence_key] = assignment
    for pool, rows in pools.items():
        signatures = {
            (
                row.get("fastq_1", ""),
                row.get("fastq_2", ""),
                row.get("sra_id", ""),
                row.get("sra_layout", ""),
            )
            for row in rows
        }
        if len(signatures) != 1:
            errors.append(f"Rows in pool {pool} must repeat identical FASTQ/SRA inputs")
        if (
            len(rows) > 1
            and stage == "raw"
            and _enabled(config, "inline_barcode_demultiplexing")
        ):
            if any(not row.get("inline_barcodes") for row in rows):
                errors.append(
                    f"Every row in shared pool {pool} requires inline_barcodes"
                )

    for module_name in ("coverage", "crosslink_sites", "enrichment"):
        if not _enabled(config, module_name):
            continue
        settings = _module(config, module_name)
        if str(settings.get("read", "any")) not in {"R1", "R2", "any"}:
            errors.append(f"modules.{module_name}.read must be R1, R2, or any")
        if str(settings.get("strand", "alignment")) not in {"alignment", "reverse"}:
            errors.append(f"modules.{module_name}.strand must be alignment or reverse")
    if _enabled(config, "crosslink_sites") and str(
        _module(config, "crosslink_sites").get("end", "5prime")
    ) not in {"5prime", "3prime"}:
        errors.append("modules.crosslink_sites.end must be 5prime or 3prime")


def validate_resolved_config(
    config: dict[str, Any], samples: list[dict[str, str]]
) -> tuple[list[str], list[str]]:
    """Return contextual validation errors and warnings for a v2 run."""
    errors: list[str] = []
    warnings: list[str] = []
    assay = config.get("assay")

    if config.get("schema_version") != 2:
        errors.append("Run config must set schema_version: 2")
    if assay not in {"generic", "rnaseq", "atacseq", "eclip"}:
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
    _validate_eclip(config, samples, errors, warnings)

    active_tools = config.get("tools", {}) or {}
    missing_tools = [
        tool for tool in _required_tools(config, samples) if tool not in active_tools
    ]
    if missing_tools:
        errors.append(
            "Missing tool configurations: "
            + ", ".join(missing_tools)
            + ". Add each with: python config/bin/ngsflow.py add-tool <tool>"
        )

    reference = config.get("reference", {}) or {}
    if (
        aligner in INDEX_KEYS
        and index_requires_build(config, aligner)
        and not genome_fasta(config)
    ):
        errors.append(
            f"reference.fasta is required because the {aligner} index must be built"
        )

    if assay == "rnaseq" and _enabled(config, "gene_counts"):
        modes = _module(config, "gene_counts").get("modes", {}) or {}
        if not any(bool(value) for value in modes.values()):
            errors.append(
                "modules.gene_counts is enabled but no counting mode is enabled"
            )
        if not reference.get("gtf"):
            errors.append("reference.gtf is required by modules.gene_counts")
        if (
            _module(config, "gene_counts").get("tool", "featurecounts")
            != "featurecounts"
        ):
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
    if _enabled(config, "rsem_quantification"):
        rsem_quant = (active_tools.get("rsem", {}).get("params", {}) or {}).get(
            "quant", {}
        ) or {}
        conflicts = sorted(RSEM_RUN_KEYS.intersection(rsem_quant))
        if conflicts:
            errors.append(
                "RSEM layout, strandedness, and execution inputs are inferred by the workflow; "
                "remove from tools/rsem.yaml params.quant: " + ", ".join(conflicts)
            )
    if assay == "rnaseq" and transcript_quant:
        if not reference.get("gtf"):
            errors.append("reference.gtf is required by transcript quantification")
        if not reference.get("fasta"):
            errors.append("reference.fasta is required by transcript quantification")

    needs_transcriptome = (
        bool(_module(config, "alignment").get("transcriptome_bam")) or transcript_quant
    )
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
        replicate_values = {
            row.get("replicate") for row in samples if row.get("replicate")
        }
        if len(replicate_values) < 2:
            warnings.append(
                "ENCODE ATAC-seq standards recommend two or more biological replicates"
            )
        if _enabled(config, "bam_filter") and not reference.get("blacklist"):
            warnings.append(
                "ATAC-seq blacklist filtering is recommended; reference.blacklist is blank"
            )
        if _enabled(config, "tss_enrichment"):
            if not reference.get("tss_bed"):
                errors.append("reference.tss_bed is required by modules.tss_enrichment")
            if not reference.get("chrom_sizes"):
                errors.append(
                    "reference.chrom_sizes is required by modules.tss_enrichment"
                )
        if _enabled(config, "frip") and not _enabled(config, "peak_calling"):
            errors.append("modules.frip requires modules.peak_calling")

    for ref_key in ("fasta", "gtf", "chrom_sizes", "blacklist", "tss_bed"):
        value = reference.get(ref_key)
        if value and not path_exists(
            config.get("_ngsflow", {}).get("project_root", "."), str(value)
        ):
            warnings.append(
                f"Reference path does not exist yet: reference.{ref_key} = {value}"
            )

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
            errors.append(
                f"Sample {sample} has invalid sra_layout '{row.get('sra_layout')}'"
            )
        for fastq_col in ("fastq_1", "fastq_2"):
            fastq = row.get(fastq_col)
            if fastq and not _project_path(config, fastq).exists():
                warnings.append(f"FASTQ path does not exist yet: {fastq}")

    return errors, warnings
