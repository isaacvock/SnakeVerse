from __future__ import annotations

import csv
import io
from collections import defaultdict
from pathlib import Path
from typing import Any


GENERIC_REQUIRED_COLUMNS = ["sample_id", "unit_id", "fastq_1"]
RNASEQ_REQUIRED_COLUMNS = [
    "sample_id",
    "unit_id",
    "fastq_1",
    "condition",
    "replicate",
    "strandedness",
]
ATACSEQ_REQUIRED_COLUMNS = [
    "sample_id",
    "unit_id",
    "fastq_1",
    "condition",
    "replicate",
]


def required_columns_for_assay(assay: str) -> list[str]:
    if assay == "rnaseq":
        return RNASEQ_REQUIRED_COLUMNS
    if assay == "atacseq":
        return ATACSEQ_REQUIRED_COLUMNS
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
        text = "".join(line for line in handle if not line.lstrip().startswith("#"))
        reader = csv.DictReader(io.StringIO(text), delimiter="\t")
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
        for row in rows:
            row.setdefault("fastq_2", "")
            row.setdefault("sra_id", "")
            row.setdefault("sra_layout", "")

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


def has_sra_samples(samples: list[dict[str, str]]) -> bool:
    return any(row.get("sra_id") and not row.get("fastq_1") for row in samples)


def sra_fastq_path(row: dict[str, str], results_dir: str, read: str) -> str:
    suffix = "R1" if read == "R1" else "R2"
    return f"{results_dir}/fastq/sra/{unit_key(row)}_{suffix}.fastq.gz"


def row_is_paired(row: dict[str, str]) -> bool:
    return bool(row.get("fastq_2")) or row.get("sra_layout") == "paired"


def unit_layout(samples: list[dict[str, str]], unit: str) -> str:
    return "paired" if row_is_paired(sample_by_unit(samples, unit)) else "single"


def sample_layout(samples: list[dict[str, str]], sample: str) -> str:
    layouts = {"paired" if row_is_paired(row) else "single" for row in units_by_sample(samples)[sample]}
    if len(layouts) != 1:
        raise ValueError(f"Sample {sample} mixes paired-end and single-end units")
    return layouts.pop()


def run_layout(samples: list[dict[str, str]]) -> str:
    layouts = {sample_layout(samples, sample) for sample in sample_ids(samples)}
    if len(layouts) != 1:
        return "mixed"
    return layouts.pop()


def fastq_for_read(
    samples: list[dict[str, str]], unit: str, read: str, results_dir: str | None = None
) -> str:
    row = sample_by_unit(samples, unit)
    column = "fastq_1" if read == "R1" else "fastq_2"
    value = row.get(column, "")
    if value:
        return value
    if row.get("sra_id"):
        if read == "R2" and not row_is_paired(row):
            raise ValueError(f"Sample unit {unit} is single-end and has no {column}")
        if results_dir is None:
            raise ValueError(
                f"Sample unit {unit} uses sra_id and needs results_dir to resolve {column}"
            )
        return sra_fastq_path(row, results_dir, read)
    if not value:
        raise ValueError(f"Sample unit {unit} is missing {column}")


def optional_fastq_for_read(
    samples: list[dict[str, str]], unit: str, read: str, results_dir: str | None = None
) -> list[str]:
    row = sample_by_unit(samples, unit)
    value = row.get("fastq_1" if read == "R1" else "fastq_2", "")
    if value:
        return [value]
    if row.get("sra_id") and (read == "R1" or row_is_paired(row)):
        if results_dir is None:
            raise ValueError(
                f"Sample unit {unit} uses sra_id and needs results_dir to resolve {read}"
            )
        return [sra_fastq_path(row, results_dir, read)]
    return [value] if value else []


def fastqc_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    targets: list[str] = []
    for row in samples:
        unit = unit_key(row)
        if row.get("fastq_1") or row.get("sra_id"):
            targets.append(f"{results_dir}/qc/fastqc/{unit}.R1")
        if row_is_paired(row):
            targets.append(f"{results_dir}/qc/fastqc/{unit}.R2")
    return targets


def trimmed_fastq_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    targets: list[str] = []
    for row in samples:
        unit = unit_key(row)
        targets.append(f"{results_dir}/fastq/trimmed/{unit}_R1.fastq.gz")
        if row_is_paired(row):
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


def narrowpeak_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    return [f"{results_dir}/peaks/macs3/{sample}/{sample}_peaks.narrowPeak" for sample in sample_ids(samples)]


def frip_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    return [f"{results_dir}/qc/atac/{sample}.frip.txt" for sample in sample_ids(samples)]


def library_complexity_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    return [f"{results_dir}/qc/atac/{sample}.library_complexity.txt" for sample in sample_ids(samples)]


def fragment_length_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    return [f"{results_dir}/qc/atac/{sample}.fragment_lengths.txt" for sample in sample_ids(samples)]


def tss_enrichment_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    return [f"{results_dir}/qc/atac/{sample}.tss_enrichment.txt" for sample in sample_ids(samples)]


def transcriptome_bam_targets(samples: list[dict[str, str]], results_dir: str) -> list[str]:
    return [f"{results_dir}/bam/transcriptome/{sample}.bam" for sample in sample_ids(samples)]


def alignment_fastq_1(
    samples: list[dict[str, str]], sample: str, results_dir: str, use_trimmed: bool
) -> list[str]:
    rows = units_by_sample(samples)[sample]
    if use_trimmed:
        return [f"{results_dir}/fastq/trimmed/{unit_key(row)}_R1.fastq.gz" for row in rows]
    return [fastq_for_read(samples, unit_key(row), "R1", results_dir) for row in rows]


def alignment_fastq_2(
    samples: list[dict[str, str]], sample: str, results_dir: str, use_trimmed: bool
) -> list[str]:
    rows = units_by_sample(samples)[sample]
    if use_trimmed:
        return [
            f"{results_dir}/fastq/trimmed/{unit_key(row)}_R2.fastq.gz"
            for row in rows
            if row_is_paired(row)
        ]
    return [
        fastq_for_read(samples, unit_key(row), "R2", results_dir)
        for row in rows
        if row_is_paired(row)
    ]


def sample_strandedness(samples: list[dict[str, str]], sample: str) -> str:
    values = {row.get("strandedness", "unstranded") for row in units_by_sample(samples)[sample]}
    if len(values) != 1:
        raise ValueError(f"Sample {sample} has mixed strandedness values: {sorted(values)}")
    return values.pop() or "unstranded"


def as_csv(values: Any) -> str:
    if isinstance(values, str):
        return values
    return ",".join(str(value) for value in values)


def as_space(values: Any) -> str:
    if isinstance(values, str):
        return values
    return " ".join(str(value) for value in values)


def bowtie2_reads_arg(r1: Any, r2: Any) -> str:
    if r2:
        return f"-1 {as_csv(r1)} -2 {as_csv(r2)}"
    return f"-U {as_csv(r1)}"


def star_reads_arg(r1: Any, r2: Any) -> str:
    if r2:
        return f"{as_csv(r1)} {as_csv(r2)}"
    return as_csv(r1)


def bwa_mem2_reads_arg(r1: Any, r2: Any) -> str:
    r1_files = as_space(r1)
    if r2:
        return f"<(gzip -cdf {r1_files}) <(gzip -cdf {as_space(r2)})"
    return f"<(gzip -cdf {r1_files})"


def featurecounts_paired_end(samples: list[dict[str, str]], config: dict[str, Any]) -> bool:
    params = config.get("tools", {}).get("featurecounts", {}).get("params", {})
    value = params.get("paired_end", "auto")
    layout = run_layout(samples)
    if layout == "mixed":
        raise ValueError("featureCounts requires all samples in a run to share PE/SE layout")
    if value == "auto":
        return layout == "paired"
    if bool(value) and layout == "single":
        raise ValueError("featureCounts paired_end is true, but the sample sheet is single-end")
    return bool(value)


def featurecounts_count_read_pairs(samples: list[dict[str, str]], config: dict[str, Any]) -> bool:
    params = config.get("tools", {}).get("featurecounts", {}).get("params", {})
    value = params.get("count_read_pairs", "auto")
    if value == "auto":
        return featurecounts_paired_end(samples, config)
    return bool(value)
