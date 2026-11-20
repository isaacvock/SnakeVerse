#!/usr/bin/env python3
"""Create deterministic tiny raw and pre-extracted eCLIP/ecCLIP FASTQs for tests."""

from __future__ import annotations

import gzip
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / ".test" / "generated" / "eclip"

ILLUMINA_R1_3P = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ILLUMINA_R2_3P = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


def write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="ascii") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def plus_strand_eclip_pair(genome: str, start: int) -> tuple[str, str]:
    """Return (R1_antisense, R2_sense) for a + strand transcript fragment."""
    r2_sense = genome[start : start + 50]
    r1_antisense = reverse_complement(genome[start + 100 : start + 150])
    return r1_antisense, r2_sense


def minus_strand_eclip_pair(genome: str, start: int) -> tuple[str, str]:
    """Return (R1_antisense, R2_sense) for a - strand transcript fragment."""
    r2_sense = reverse_complement(genome[start + 100 : start + 150])
    r1_antisense = genome[start : start + 50]
    return r1_antisense, r2_sense


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
    plus_0 = plus_strand_eclip_pair(genome, 5_000)
    plus_1 = plus_strand_eclip_pair(genome, 10_000)
    minus_0 = minus_strand_eclip_pair(genome, 15_000)

    # A01 = AAGCAAT (revcomp ATTGCTT), B06 = GGCTTGT (revcomp ACAAGCC).
    # Includes within-A01 PCR duplicate, cross-barcode identical-locus + identical-UMI
    # collision (must survive --per-cell --cell-tag=BC), 1-mismatch UMI (must survive
    # method: unique), and double-ligation artifacts on both R1 and R2.
    encode = [
        ("dup_a01_1", "AAGCAAT" + plus_0[0], "AAAAA" + plus_0[1]),
        ("dup_a01_2", "AAGCAAT" + plus_0[0], "AAAAA" + plus_0[1]),
        ("same_umi_b06", "GGCTTGT" + plus_0[0], "AAAAA" + plus_0[1]),
        ("umi_1mm_a01", "AAGCAAT" + plus_0[0], "AAAAC" + plus_0[1]),
        (
            "dbl_lig_b06",
            "GGCTTGT" + "CTTCCGATCTGGCTTGT" + minus_0[0] + "GGGGG" + ILLUMINA_R1_3P[:20],
            "CCCCC"
            + minus_0[1]
            + "CAAGCC"
            + ILLUMINA_R2_3P[:18]
            + "ACAAGCC"
            + ILLUMINA_R2_3P[:15],
        ),
    ]
    write_fastq(OUT / "encode_pool_R1.fastq.gz", [(name, r1) for name, r1, _ in encode])
    write_fastq(OUT / "encode_pool_R2.fastq.gz", [(name, r2) for name, _, r2 in encode])

    smi = [
        ("smi_plus", plus_1[0], "TTTTT" + plus_1[1]),
        ("smi_minus", minus_0[0], "GGGGG" + minus_0[1]),
    ]
    write_fastq(OUT / "encode_smi_R1.fastq.gz", [(name, r1) for name, r1, _ in smi])
    write_fastq(OUT / "encode_smi_R2.fastq.gz", [(name, r2) for name, _, r2 in smi])

    # ENCODE DCC pre-extracted style (<UMI>:<read_name>) split across two barcoded units
    # sharing identical coordinates + UMI AAAAA.
    pre_a01 = [
        ("AAAAA:pre_a01_1", plus_0[0], plus_0[1]),
        ("AAAAA:pre_a01_2", plus_0[0], plus_0[1]),
    ]
    pre_b06 = [
        ("AAAAA:pre_b06_1", plus_0[0], plus_0[1]),
    ]
    write_fastq(OUT / "preextracted_A01_R1.fastq.gz", [(n, r1) for n, r1, _ in pre_a01])
    write_fastq(OUT / "preextracted_A01_R2.fastq.gz", [(n, r2) for n, _, r2 in pre_a01])
    write_fastq(OUT / "preextracted_B06_R1.fastq.gz", [(n, r1) for n, r1, _ in pre_b06])
    write_fastq(OUT / "preextracted_B06_R2.fastq.gz", [(n, r2) for n, _, r2 in pre_b06])

    # ecCLIP: R1 sense (+3 nt 5' overhang + noisy poly(A) tail), R2 8-nt UMI-only.
    ec_plus_r1 = genome[5_000 : 5_050]
    ec_minus_r1 = reverse_complement(genome[15_000 : 15_050])
    noisy_polya = "A" * 10 + "C" + "A" * 12
    ecclip = [
        ("ec_dup_a", "GGA" + ec_plus_r1 + "A" * 18, "ACGTACGT"),
        ("ec_dup_b", "GGA" + ec_plus_r1 + "A" * 18, "ACGTACGT"),
        ("ec_umi_b", "GGA" + ec_plus_r1 + noisy_polya, "TGCATGCA"),
        ("ec_minus_1", "GGA" + ec_minus_r1 + "A" * 18, "GGGGCCCC"),
    ]
    write_fastq(OUT / "ecclip_ip_R1.fastq.gz", [(name, r1) for name, r1, _ in ecclip])
    write_fastq(OUT / "ecclip_ip_R2.fastq.gz", [(name, r2) for name, _, r2 in ecclip])

    tag = [("tag_1", "GGA" + genome[10_000 : 10_050] + "A" * 18, "AAAACCCC")]
    write_fastq(OUT / "ecclip_tag_R1.fastq.gz", [(name, r1) for name, r1, _ in tag])
    write_fastq(OUT / "ecclip_tag_R2.fastq.gz", [(name, r2) for name, _, r2 in tag])


if __name__ == "__main__":
    main()
