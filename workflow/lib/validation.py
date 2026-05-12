from __future__ import annotations

from pathlib import Path
from typing import Any

from samples import required_columns_for_assay


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
    if assay not in {"generic", "rnaseq"}:
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
    if config.get("steps", {}).get("coverage", False):
        required_tools.append("deeptools")
    for tool in sorted(set(required_tools)):
        if tool not in config.get("tools", {}):
            missing_tools.append(tool)
    if missing_tools:
        errors.append("Missing active tool profiles: " + ", ".join(missing_tools))

    genome = config.get("genome", {}) or {}
    if aligner == "bowtie2" and not genome.get("bowtie2_index"):
        errors.append("Genome profile must define genome.bowtie2_index for Bowtie2 runs")
    if aligner == "star" and not genome.get("star_index"):
        errors.append("Genome profile must define genome.star_index for STAR runs")
    if assay == "rnaseq" and not genome.get("gtf"):
        errors.append("Genome profile must define genome.gtf for RNA-seq featureCounts")

    for ref_key in ("fasta", "gtf", "chrom_sizes", "bowtie2_index", "star_index"):
        value = genome.get(ref_key)
        if value:
            ref_path = _project_path(config, str(value))
            if not ref_path.exists():
                warnings.append(f"Reference path for genome.{ref_key} does not exist yet: {value}")

    required_cols = required_columns_for_assay(str(assay))
    for row in samples:
        for column in required_cols:
            if not row.get(column):
                errors.append(f"Sample {row.get('sample_id', '<unknown>')} is missing {column}")
        for fastq_col in ("fastq_1", "fastq_2"):
            fastq = row.get(fastq_col)
            if fastq and not _project_path(config, fastq).exists():
                warnings.append(f"FASTQ path does not exist yet: {fastq}")

    return errors, warnings

