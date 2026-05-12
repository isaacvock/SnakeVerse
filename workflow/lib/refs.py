from __future__ import annotations

from pathlib import Path
from typing import Any


INDEX_KEYS = {
    "bowtie2": "bowtie2_index",
    "star": "star_index",
    "bwa_mem2": "bwa_mem2_index",
}


def genome_name(config: dict[str, Any]) -> str:
    return str(config.get("genome", {}).get("name") or "genome")


def genome_slug(config: dict[str, Any]) -> str:
    return "".join(char if char.isalnum() or char in ("-", "_", ".") else "_" for char in genome_name(config))


def configured_index(config: dict[str, Any], aligner: str) -> str:
    key = INDEX_KEYS.get(aligner)
    if not key:
        return ""
    return str(config.get("genome", {}).get(key) or "").rstrip("/")


def index_is_configured(config: dict[str, Any], aligner: str) -> bool:
    return bool(configured_index(config, aligner))


def generated_index_dir(config: dict[str, Any], results_dir: str, aligner: str) -> str:
    return f"{results_dir}/reference/{aligner}/{genome_slug(config)}"


def generated_index_prefix(config: dict[str, Any], results_dir: str, aligner: str) -> str:
    index_dir = generated_index_dir(config, results_dir, aligner)
    if aligner in {"bowtie2", "bwa_mem2"}:
        return f"{index_dir}/{genome_slug(config)}"
    return index_dir


def aligner_index_prefix(config: dict[str, Any], results_dir: str, aligner: str) -> str:
    return configured_index(config, aligner) or generated_index_prefix(config, results_dir, aligner)


def generated_index_marker(config: dict[str, Any], results_dir: str, aligner: str) -> str:
    return f"{generated_index_dir(config, results_dir, aligner)}/.snakeverse_{aligner}_index.done"


def aligner_index_inputs(config: dict[str, Any], results_dir: str, aligner: str) -> list[str]:
    if index_is_configured(config, aligner):
        return []
    return [generated_index_marker(config, results_dir, aligner)]


def genome_fasta(config: dict[str, Any]) -> str:
    return str(config.get("genome", {}).get("fasta") or "")


def genome_gtf(config: dict[str, Any]) -> str:
    return str(config.get("genome", {}).get("gtf") or "")


def star_gtf_arg(config: dict[str, Any]) -> str:
    gtf = genome_gtf(config)
    return f"--sjdbGTFfile {gtf}" if gtf else ""


def path_exists(project_root: str | Path, value: str) -> bool:
    path = Path(value)
    if not path.is_absolute():
        path = Path(project_root) / path
    return path.exists()
