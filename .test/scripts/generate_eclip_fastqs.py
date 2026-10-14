#!/usr/bin/env python3
"""Create deterministic tiny raw eCLIP/ecCLIP FASTQs for integration tests."""

from __future__ import annotations

import gzip
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / ".test" / "generated" / "eclip"


def write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="ascii") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_reference() -> str:
    """Write a compact, non-palindromic reference with unique test loci."""
    genome = "".join(random.Random(20261014).choices("ACGT", k=200_000))
    fasta = OUT / "genome.fasta"
    fasta.parent.mkdir(parents=True, exist_ok=True)
    with fasta.open("w", encoding="ascii") as handle:
        handle.write(">chrClip\n")
        for offset in range(0, len(genome), 60):
            handle.write(genome[offset : offset + 60] + "\n")
    (OUT / "genome.chrom.sizes").write_text(
        f"chrClip\t{len(genome)}\n", encoding="ascii"
    )
    (OUT / "genome.gtf").write_text(
        'chrClip\tfixture\texon\t1001\t2000\t.\t+\t.\tgene_id "clip_gene"; '
        'transcript_id "clip_tx";\n',
        encoding="ascii",
    )
    return genome


def main() -> None:
    genome = write_reference()
    loci = [
        (genome[start : start + 50], reverse_complement(genome[start + 100 : start + 150]))
        for start in (5_000, 10_000, 15_000)
    ]

    encode = [
        ("dup_a", "ATTGCTT" + loci[0][0], "AAAAA" + loci[0][1]),
        ("dup_b", "ATTGCTT" + loci[0][0], "AAAAA" + loci[0][1]),
        ("umi_b", "ATTGCTT" + loci[0][0], "CCCCC" + loci[0][1]),
        ("barcode_b", "ACAAGCC" + loci[1][0], "GGGGG" + loci[1][1]),
    ]
    write_fastq(OUT / "encode_pool_R1.fastq.gz", [(name, r1) for name, r1, _ in encode])
    write_fastq(OUT / "encode_pool_R2.fastq.gz", [(name, r2) for name, _, r2 in encode])

    smi = [("smi_1", loci[2][0], "TTTTT" + loci[2][1])]
    write_fastq(OUT / "encode_smi_R1.fastq.gz", [(name, r1) for name, r1, _ in smi])
    write_fastq(OUT / "encode_smi_R2.fastq.gz", [(name, r2) for name, _, r2 in smi])

    preextracted = [
        ("pre_dup_a|UMI=AAAAA", loci[0][0], loci[0][1]),
        ("pre_dup_b|UMI=AAAAA", loci[0][0], loci[0][1]),
        ("pre_umi_b|UMI=CCCCC", loci[0][0], loci[0][1]),
    ]
    write_fastq(
        OUT / "preextracted_R1.fastq.gz",
        [(name, r1) for name, r1, _ in preextracted],
    )
    write_fastq(
        OUT / "preextracted_R2.fastq.gz",
        [(name, r2) for name, _, r2 in preextracted],
    )

    ecclip = [
        ("ec_dup_a", "GGA" + loci[0][0] + "A" * 18, "ACGTACGT"),
        ("ec_dup_b", "GGA" + loci[0][0] + "A" * 18, "ACGTACGT"),
        ("ec_umi_b", "GGA" + loci[0][0] + "A" * 18, "TGCATGCA"),
    ]
    write_fastq(OUT / "ecclip_ip_R1.fastq.gz", [(name, r1) for name, r1, _ in ecclip])
    write_fastq(OUT / "ecclip_ip_R2.fastq.gz", [(name, r2) for name, _, r2 in ecclip])

    tag = [("tag_1", "GGA" + loci[1][0] + "A" * 18, "AAAACCCC")]
    write_fastq(OUT / "ecclip_tag_R1.fastq.gz", [(name, r1) for name, r1, _ in tag])
    write_fastq(OUT / "ecclip_tag_R2.fastq.gz", [(name, r2) for name, _, r2 in tag])


if __name__ == "__main__":
    main()
