from __future__ import annotations

from pathlib import Path
from typing import Any


INDEX_KEYS = {
    "bowtie2": "bowtie2",
    "star": "star",
    "bwa_mem2": "bwa_mem2",
}

STAR_INDEX_FILES = (
    "Genome",
    "SA",
    "SAindex",
    "chrLength.txt",
    "chrName.txt",
    "chrNameLength.txt",
    "chrStart.txt",
    "genomeParameters.txt",
)


def genome_name(config: dict[str, Any]) -> str:
    return str(config.get("reference", {}).get("name") or "genome")


def genome_slug(config: dict[str, Any]) -> str:
    return "".join(char if char.isalnum() or char in ("-", "_", ".") else "_" for char in genome_name(config))


def configured_index(config: dict[str, Any], aligner: str) -> str:
    key = INDEX_KEYS.get(aligner)
    if not key:
        return ""
    indexes = config.get("reference", {}).get("indexes", {}) or {}
    return str(indexes.get(key) or "").rstrip("/")


def index_is_configured(config: dict[str, Any], aligner: str) -> bool:
    return bool(configured_index(config, aligner))


def _project_path(config: dict[str, Any], value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    project_root = Path(config.get("_ngsflow", {}).get("project_root", "."))
    return project_root / path


def configured_index_is_ready(config: dict[str, Any], aligner: str) -> bool:
    """Return true when a configured index contains the expected files."""
    configured = configured_index(config, aligner)
    if not configured:
        return False

    index_path = _project_path(config, configured)
    if aligner == "star":
        return all((index_path / name).exists() for name in STAR_INDEX_FILES)
    if aligner == "bowtie2":
        small_index = (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
        large_index = tuple(suffix + "l" for suffix in small_index)
        return all(Path(str(index_path) + suffix).exists() for suffix in small_index) or all(
            Path(str(index_path) + suffix).exists() for suffix in large_index
        )
    if aligner == "bwa_mem2":
        suffixes = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
        return all(Path(str(index_path) + suffix).exists() for suffix in suffixes)
    return index_path.exists()


def index_requires_build(config: dict[str, Any], aligner: str) -> bool:
    return not configured_index_is_ready(config, aligner)


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


def aligner_index_marker(config: dict[str, Any], results_dir: str, aligner: str) -> str:
    configured = configured_index(config, aligner)
    if not configured:
        return generated_index_marker(config, results_dir, aligner)
    if aligner == "star":
        return f"{configured}/.snakeverse_star_index.done"
    return f"{configured}.snakeverse_{aligner}_index.done"


def aligner_index_inputs(config: dict[str, Any], results_dir: str, aligner: str) -> list[str]:
    if configured_index_is_ready(config, aligner):
        return []
    return [aligner_index_marker(config, results_dir, aligner)]


def genome_fasta(config: dict[str, Any]) -> str:
    return str(config.get("reference", {}).get("fasta") or "")


def genome_gtf(config: dict[str, Any]) -> str:
    return str(config.get("reference", {}).get("gtf") or "")


def star_gtf_arg(config: dict[str, Any]) -> str:
    gtf = genome_gtf(config)
    return f"--sjdbGTFfile {gtf}" if gtf else ""


def path_exists(project_root: str | Path, value: str) -> bool:
    path = Path(value)
    if not path.is_absolute():
        path = Path(project_root) / path
    return path.exists()
