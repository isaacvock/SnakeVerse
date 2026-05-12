from __future__ import annotations

from pathlib import Path
from typing import Any

from refs import INDEX_KEYS, configured_index, genome_fasta, path_exists
from samples import featurecounts_paired_end, required_columns_for_assay, sample_ids, sample_layout


def _project_path(config: dict[str, Any], value: str) -> Path:
    root = Path(config.get("_ngsflow", {}).get("project_root", "."))
    path = Path(value)
    if path.is_absolute():
        return path
    return root / path


def validate_resolved_config(config: dict[str, Any], samples: list[dict[str, str]]) -> tuple[list[str], list[str]]:
    """Return validation errors and warnings for a resolved workflow config."""
    errors: list[str] = []
    warnings: list[str] = []

    assay = config.get("assay")
    if assay not in {"generic", "rnaseq", "atacseq"}:
        errors.append(f"Unsupported assay: {assay}")

    for key in ("project", "samples", "results_dir"):
        if key not in config:
            errors.append(f"Resolved config is missing required key: {key}")

    missing_tools = []
    aligner = config.get("alignment", {}).get("tool")
    required_tools = ["fastqc", "samtools", "multiqc"]
    if config.get("steps", {}).get("trimming", False):
        required_tools.append("cutadapt")
    if aligner:
        required_tools.append(aligner)
    if assay == "rnaseq":
        required_tools.append("featurecounts")
    if assay == "atacseq":
        required_tools.extend(["macs3", "bedtools"])
    if config.get("steps", {}).get("coverage", False):
        required_tools.append("deeptools")
    for tool in sorted(set(required_tools)):
        if tool not in config.get("tools", {}):
            missing_tools.append(tool)
    if missing_tools:
        errors.append("Missing active tool profiles: " + ", ".join(missing_tools))

    genome = config.get("genome", {}) or {}
    if aligner not in {"bowtie2", "star", "bwa_mem2"}:
        errors.append(f"Unsupported alignment.tool: {aligner}")
    if aligner in INDEX_KEYS and not configured_index(config, aligner) and not genome_fasta(config):
        errors.append(f"Genome profile must define genome.fasta when building a {aligner} index")
    if assay == "rnaseq" and not genome.get("gtf"):
        errors.append("Genome profile must define genome.gtf for RNA-seq featureCounts")
    if assay == "rnaseq" and config.get("outputs", {}).get("gene_counts", False):
        try:
            featurecounts_paired_end(samples, config)
        except ValueError as exc:
            errors.append(str(exc))
    if config.get("outputs", {}).get("transcriptome_bam", False):
        if aligner != "star":
            errors.append("outputs.transcriptome_bam requires alignment.tool: star")
        star_params = config.get("tools", {}).get("star", {}).get("params", {})
        align_params = star_params.get("align", star_params)
        quant_mode = str(align_params.get("quantMode", ""))
        if "TranscriptomeSAM" not in quant_mode:
            errors.append("STAR transcriptome BAM output requires star.params.align.quantMode to include TranscriptomeSAM")
    if assay == "atacseq":
        replicate_values = {row.get("replicate") for row in samples if row.get("replicate")}
        if len(replicate_values) < 2:
            warnings.append("ENCODE ATAC-seq standards recommend two or more biological replicates")
        if not genome.get("blacklist"):
            warnings.append("ATAC-seq blacklist filtering is recommended; genome.blacklist is blank")
        if config.get("outputs", {}).get("tss_enrichment", False) and not genome.get("tss_bed"):
            errors.append("outputs.tss_enrichment requires genome.tss_bed")

    for ref_key in ("fasta", "gtf", "chrom_sizes", "blacklist", "tss_bed", "bowtie2_index", "star_index", "bwa_mem2_index"):
        value = genome.get(ref_key)
        if value:
            if not path_exists(config.get("_ngsflow", {}).get("project_root", "."), str(value)):
                warnings.append(f"Reference path for genome.{ref_key} does not exist yet: {value}")

    required_cols = required_columns_for_assay(str(assay))
    for sample in sample_ids(samples):
        try:
            sample_layout(samples, sample)
        except ValueError as exc:
            errors.append(str(exc))
    for row in samples:
        for column in required_cols:
            if not row.get(column):
                errors.append(f"Sample {row.get('sample_id', '<unknown>')} is missing {column}")
        for fastq_col in ("fastq_1", "fastq_2"):
            fastq = row.get(fastq_col)
            if fastq and not _project_path(config, fastq).exists():
                warnings.append(f"FASTQ path does not exist yet: {fastq}")

    return errors, warnings
