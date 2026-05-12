from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Any


GENERIC_REQUIRED_COLUMNS = ["sample_id", "unit_id", "fastq_1", "fastq_2"]
RNASEQ_REQUIRED_COLUMNS = [
    "sample_id",
    "unit_id",
    "fastq_1",
    "fastq_2",
    "condition",
    "replicate",
    "strandedness",
]


def required_columns_for_assay(assay: str) -> list[str]:
    if assay == "rnaseq":
        return RNASEQ_REQUIRED_COLUMNS
    return GENERIC_REQUIRED_COLUMNS


def project_path(path: str | Path, project_root: str | Path | None = None) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    root = Path(project_root) if project_root else Path.cwd()
    return (root / candidate).resolve()


def load_samples(
    path: str | Path, assay: str, project_root: str | Path | None = None
) -> list[dict[str, str]]:
    sample_path = project_path(path, project_root)
    if not sample_path.exists():
        raise FileNotFoundError(f"Samples file does not exist: {sample_path}")

    with sample_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing = [col for col in required_columns_for_assay(assay) if col not in fieldnames]
        if missing:
            raise ValueError(
                f"Samples file {sample_path} is missing required columns: "
                + ", ".join(missing)
            )
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
            if any((value or "").strip() for value in row.values())
        ]

    seen_units: set[str] = set()
    for row in rows:
        unit = unit_key(row)
        if unit in seen_units:
            raise ValueError(f"Duplicate sample/unit combination in samples file: {unit}")
        seen_units.add(unit)
    return rows


def sample_ids(samples: list[dict[str, str]]) -> list[str]:
    return sorted({row["sample_id"] for row in samples})


def unit_key(row: dict[str, str]) -> str:
    return f"{row['sample_id']}__{row['unit_id']}"


def unit_ids(samples: list[dict[str, str]]) -> list[str]:
    return [unit_key(row) for row in samples]


def units_by_sample(samples: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in samples:
        grouped[row["sample_id"]].append(row)
    return dict(grouped)


def sample_by_unit(samples: list[dict[str, str]], unit: str) -> dict[str, str]:
    for row in samples:
        if unit_key(row) == unit:
            return row
    raise KeyError(f"Unknown sample unit: {unit}")


def fastq_for_read(samples: list[dict[str, str]], unit: str, read: str) -> str:
    row = sample_by_unit(samples, unit)
    column = "fastq_1" if read == "R1" else "fastq_2"
    value = row.get(column, "")
    if not value:
        raise ValueError(f"Sample unit {unit} is missing {column}")
    return value


def fastqc_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    targets: list[str] = []
    for row in samples:
        unit = unit_key(row)
        if row.get("fastq_1"):
            targets.append(f"{results_dir}/qc/fastqc/{unit}.R1")
        if row.get("fastq_2"):
            targets.append(f"{results_dir}/qc/fastqc/{unit}.R2")
    return targets


def trimmed_fastq_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    targets: list[str] = []
    for row in samples:
        unit = unit_key(row)
        targets.append(f"{results_dir}/fastq/trimmed/{unit}_R1.fastq.gz")
        if row.get("fastq_2"):
            targets.append(f"{results_dir}/fastq/trimmed/{unit}_R2.fastq.gz")
    return targets


def raw_bam_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    return [
        item
        for sample in sample_ids(samples)
        for item in [
            f"{results_dir}/bam/raw/{sample}.bam",
            f"{results_dir}/bam/raw/{sample}.bam.bai",
        ]
    ]


def filtered_bam_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    return [
        item
        for sample in sample_ids(samples)
        for item in [
            f"{results_dir}/bam/filtered/{sample}.bam",
            f"{results_dir}/bam/filtered/{sample}.bam.bai",
        ]
    ]


def bam_qc_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    return [
        item
        for sample in sample_ids(samples)
        for item in [
            f"{results_dir}/qc/bam/{sample}.flagstat.txt",
            f"{results_dir}/qc/bam/{sample}.idxstats.txt",
        ]
    ]


def bigwig_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    return [f"{results_dir}/tracks/bigwig/{sample}.bw" for sample in sample_ids(samples)]


def alignment_fastq_1(
    samples: list[dict[str, str]], sample: str, results_dir: str, use_trimmed: bool
) -> list[str]:
    rows = units_by_sample(samples)[sample]
    if use_trimmed:
        return [f"{results_dir}/fastq/trimmed/{unit_key(row)}_R1.fastq.gz" for row in rows]
    return [row["fastq_1"] for row in rows]


def alignment_fastq_2(
    samples: list[dict[str, str]], sample: str, results_dir: str, use_trimmed: bool
) -> list[str]:
    rows = units_by_sample(samples)[sample]
    if use_trimmed:
        return [f"{results_dir}/fastq/trimmed/{unit_key(row)}_R2.fastq.gz" for row in rows]
    return [row["fastq_2"] for row in rows if row.get("fastq_2")]


def sample_strandedness(samples: list[dict[str, str]], sample: str) -> str:
    values = {row.get("strandedness", "unstranded") for row in units_by_sample(samples)[sample]}
    if len(values) != 1:
        raise ValueError(f"Sample {sample} has mixed strandedness values: {sorted(values)}")
    return values.pop() or "unstranded"


def as_csv(values: Any) -> str:
    if isinstance(values, str):
        return values
    return ",".join(str(value) for value in values)

