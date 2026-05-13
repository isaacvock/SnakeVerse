#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Span:
    chrom: str
    start: int
    end: int
    strand: str


def parse_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for field in raw.rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field and " " not in field.split("=", 1)[0]:
            key, value = field.split("=", 1)
        else:
            parts = field.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def update_span(spans: dict[str, Span], gene_id: str, chrom: str, start: int, end: int, strand: str) -> None:
    if gene_id not in spans:
        spans[gene_id] = Span(chrom=chrom, start=start, end=end, strand=strand)
        return
    span = spans[gene_id]
    if chrom != span.chrom:
        return
    span.start = min(span.start, start)
    span.end = max(span.end, end)
    if strand != span.strand:
        span.strand = "."


def build_gene_saf(gtf: Path, output: Path) -> None:
    gene_spans: dict[str, Span] = {}
    derived_spans: dict[str, Span] = {}

    with gtf.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start_text, end_text, _, strand, _, raw_attrs = fields
            attrs = parse_attributes(raw_attrs)
            gene_id = attrs.get("gene_id") or attrs.get("ID") or attrs.get("gene")
            if not gene_id:
                continue
            start = int(start_text)
            end = int(end_text)
            update_span(derived_spans, gene_id, chrom, start, end, strand)
            if feature == "gene":
                update_span(gene_spans, gene_id, chrom, start, end, strand)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        handle.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for gene_id in sorted(derived_spans, key=lambda key: (derived_spans[key].chrom, derived_spans[key].start, key)):
            span = gene_spans.get(gene_id, derived_spans[gene_id])
            handle.write(f"{gene_id}\t{span.chrom}\t{span.start}\t{span.end}\t{span.strand}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert a GTF to featureCounts SAF gene spans.")
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    build_gene_saf(args.gtf, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
