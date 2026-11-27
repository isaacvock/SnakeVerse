#!/usr/bin/env python3
"""Create deterministic synthetic FASTQ and reference fixtures for mutational_profiling."""

from __future__ import annotations

import gzip
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / ".test" / "generated" / "mutprof"


def write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="ascii") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def build_genome() -> str:
    rng = random.Random(20261126)
    bases = [rng.choice("ACGT") for _ in range(200_000)]

    # Ensure deterministic reference bases at target test sites:
    # 1) Plus-strand exon (1001..2000):
    #    - 0-based 1150 (pos 1151): ref 'A' (true RBP AG editing site in RBP_1)
    #    - 0-based 1155 (pos 1156): ref 'A' (germline A->G SNP in CTRL_1 + CTRL_2 + RBP_1)
    bases[1150] = "A"
    bases[1155] = "A"

    # 2) Minus-strand exon (5001..6000):
    #    - 0-based 5150 (pos 5151): + ref 'G' => minus-strand RNA ref 'C' (true CT edit in RBP_1)
    bases[5150] = "G"

    # 3) Unstranded locus (10001..11000):
    #    - 0-based 10130 (pos 10131): ref 'A' (50% AG edit, preserved by target_min_freq=0.90)
    #    - 0-based 10135 (pos 10136): ref 'A' (100% homozygous A->G SNP, masked)
    #    - 0-based 10140 (pos 10141): ref 'C' (50% heterozygous C->A SNP, masked)
    bases[10130] = "A"
    bases[10135] = "A"
    bases[10140] = "C"

    genome = "".join(bases)
    fasta = OUT / "genome.fasta"
    fasta.parent.mkdir(parents=True, exist_ok=True)
    with fasta.open("w", encoding="ascii") as handle:
        handle.write(">chrMut\n")
        for offset in range(0, len(genome), 60):
            handle.write(genome[offset : offset + 60] + "\n")

    (OUT / "genome.chrom.sizes").write_text(
        f"chrMut\t{len(genome)}\n", encoding="ascii"
    )
    (OUT / "genome.gtf").write_text(
        'chrMut\tfixture\texon\t1001\t2000\t.\t+\t.\tgene_id "gene_plus"; transcript_id "tx_plus";\n'
        'chrMut\tfixture\texon\t5001\t6000\t.\t-\t.\tgene_id "gene_minus"; transcript_id "tx_minus";\n'
        'chrMut\tfixture\texon\t10001\t11000\t.\t+\t.\tgene_id "gene_unstr"; transcript_id "tx_unstr";\n',
        encoding="ascii",
    )
    return genome


def reverse_stranded_plus_pair(
    genome_mut: str, start: int, frag_len: int = 100, read_len: int = 60
) -> tuple[str, str]:
    """Return (R1_antisense, R2_sense) for a '+' strand transcript in a reverse-stranded PE library."""
    r2 = genome_mut[start : start + read_len]
    r1 = reverse_complement(genome_mut[start + frag_len - read_len : start + frag_len])
    return r1, r2


def reverse_stranded_minus_pair(
    genome_mut: str, start: int, frag_len: int = 100, read_len: int = 60
) -> tuple[str, str]:
    """Return (R1_antisense_to_minus_RNA, R2_sense_to_minus_RNA) for a '-' strand transcript."""
    r1 = genome_mut[start : start + read_len]
    r2 = reverse_complement(genome_mut[start + frag_len - read_len : start + frag_len])
    return r1, r2


def mutate_genome(genome: str, edits: dict[int, str]) -> str:
    buf = list(genome)
    for pos_0, base in edits.items():
        buf[pos_0] = base
    return "".join(buf)


def main() -> None:
    genome = build_genome()

    # --- Fixture 1: Stranded PE dual-deaminase (AG + CT) with pooled negative controls ---
    rbp_pairs: list[tuple[str, str, str]] = []
    # 8 staggered unique fragments covering 1150 & 1155 in mate overlap, plus 2 PCR duplicates of k=0
    for k in range(8):
        edits = {1155: "G"} if k < 4 else {}  # germline SNP also present at 50% in RBP_1
        if k < 4:
            edits[1150] = "G"  # 4/8 = 50% true AG editing at pos 1151
        g_mod = mutate_genome(genome, edits)
        r1, r2 = reverse_stranded_plus_pair(g_mod, 1100 + k)
        rbp_pairs.append((f"rbp_plus_{k}", r1, r2))
        if k == 0:
            rbp_pairs.append(("rbp_plus_0_pcrdup1", r1, r2))
            rbp_pairs.append(("rbp_plus_0_pcrdup2", r1, r2))

    # 6 staggered unique fragments on '-' strand exon covering 5150 (pos 5151)
    for k in range(6):
        # '-' RNA C->T edit corresponds to + genomic G->A at 0-based 5150
        edits = {5150: "A"} if k < 3 else {}
        g_mod = mutate_genome(genome, edits)
        r1, r2 = reverse_stranded_minus_pair(g_mod, 5100 + k)
        rbp_pairs.append((f"rbp_minus_{k}", r1, r2))

    write_fastq(OUT / "rbp1_R1.fastq.gz", [(n, r1) for n, r1, _ in rbp_pairs])
    write_fastq(OUT / "rbp1_R2.fastq.gz", [(n, r2) for n, _, r2 in rbp_pairs])

    # Two negative controls: each has 2/4 reads with A->G at 1155 (pos 1156); pooled = 4/8 -> called SNP!
    for ctrl_idx in (1, 2):
        ctrl_pairs: list[tuple[str, str, str]] = []
        offset_base = 1100 if ctrl_idx == 1 else 1104
        for k in range(4):
            edits = {1155: "G"} if k < 2 else {}
            g_mod = mutate_genome(genome, edits)
            r1, r2 = reverse_stranded_plus_pair(g_mod, offset_base + k)
            ctrl_pairs.append((f"ctrl{ctrl_idx}_{k}", r1, r2))
        write_fastq(OUT / f"ctrl{ctrl_idx}_R1.fastq.gz", [(n, r1) for n, r1, _ in ctrl_pairs])
        write_fastq(OUT / f"ctrl{ctrl_idx}_R2.fastq.gz", [(n, r2) for n, _, r2 in ctrl_pairs])

    # --- Fixture 2: Unstranded SE HyperTRIBE (AG) with NO negative controls ---
    se_reads: list[tuple[str, str]] = []
    for k in range(10):
        edits = {10135: "G"}  # 10/10 = 100% homozygous AG SNP at pos 10136 (masked >= 0.90)
        if k < 5:
            edits[10130] = "G"  # 5/10 = 50% true AG edit at pos 10131 (kept < 0.90)
            edits[10140] = "A"  # 5/10 = 50% non-target C->A SNP at pos 10141 (masked >= 0.25)
        g_mod = mutate_genome(genome, edits)
        seq = g_mod[10100 + k : 10100 + k + 70]
        se_reads.append((f"unstr_{k}", seq))

    write_fastq(OUT / "unstr_rbp_R1.fastq.gz", se_reads)


if __name__ == "__main__":
    main()
