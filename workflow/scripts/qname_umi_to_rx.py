#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pysam

BARCODE_STAMP_RE = re.compile(
    r"\|BC=(?P<label>[A-Za-z0-9_.-]+):(?P<seq>[ACGTRYSWKMBDHVN]+)(?=_|$)"
)
SAM_TAG_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]$")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Copy a QNAME-encoded UMI into a SAM tag (default RX) and, when "
            "present, copy stamped inline barcode sequences into SAM tag BC."
        )
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--regex", default=r"_(?P<umi>[ACGTN]+)$")
    parser.add_argument("--tag", default="RX")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not SAM_TAG_RE.fullmatch(args.tag):
        raise SystemExit(f"--tag must be a two-character SAM tag, got: {args.tag}")
    pattern = re.compile(args.regex)
    if "umi" not in pattern.groupindex:
        raise SystemExit("--regex must contain a named group '(?P<umi>...)'")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    tagged = 0
    barcoded = 0
    with pysam.AlignmentFile(args.input, "rb") as source:
        with pysam.AlignmentFile(args.output, "wb", template=source) as destination:
            for record in source.fetch(until_eof=True):
                qname = record.query_name
                bc_match = BARCODE_STAMP_RE.search(qname)
                if bc_match is not None:
                    bc_label = bc_match.group("label")
                    bc_seq = bc_match.group("seq")
                    clean_qname = qname[: bc_match.start()] + qname[bc_match.end() :]
                    retained_qname = (
                        qname[: bc_match.start()]
                        + f"|BC={bc_label}"
                        + qname[bc_match.end() :]
                    )
                else:
                    bc_seq = None
                    clean_qname = qname
                    retained_qname = qname

                match = pattern.search(clean_qname)
                if match is None:
                    raise SystemExit(
                        f"Read name does not match UMI regex: {record.query_name}"
                    )
                if retained_qname != qname:
                    record.query_name = retained_qname
                record.set_tag(args.tag, match.group("umi"), value_type="Z", replace=True)
                if bc_seq is not None:
                    record.set_tag("BC", bc_seq, value_type="Z", replace=True)
                    barcoded += 1
                destination.write(record)
                tagged += 1
    print(
        f"Tagged {tagged} BAM records with {args.tag}"
        + (f" and {barcoded} with BC" if barcoded else "")
    )


if __name__ == "__main__":
    main()
