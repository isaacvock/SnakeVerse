from __future__ import annotations

from typing import Any

from samples import sample_strandedness


def coverage_bam_for_sample(config: dict[str, Any], results_dir: str, sample: str) -> str:
    if config.get("steps", {}).get("bam_filter", True):
        return f"{results_dir}/bam/filtered/{sample}.bam"
    return f"{results_dir}/bam/raw/{sample}.bam"


def coverage_strand_arg(config: dict[str, Any], samples: list[dict[str, str]], sample: str) -> str:
    coverage = config.get("coverage", {}) or {}
    mode = coverage.get("rna_strand_mode", "none")
    if config.get("assay") != "rnaseq" or mode in ("none", None, False):
        return ""
    if mode == "sample":
        strandedness = sample_strandedness(samples, sample)
        if strandedness in ("forward", "reverse"):
            return f"--filterRNAstrand {strandedness}"
        return ""
    if mode in ("forward", "reverse"):
        return f"--filterRNAstrand {mode}"
    return ""

