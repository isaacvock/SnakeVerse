#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pysam


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Copy a QNAME-encoded UMI into the BAM RX tag."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--regex", default=r"_(?P<umi>[ACGTN]+)$")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    pattern = re.compile(args.regex)
    if "umi" not in pattern.groupindex:
        raise SystemExit("--regex must contain a named group '(?P<umi>...)'")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    tagged = 0
    with pysam.AlignmentFile(args.input, "rb") as source:
        with pysam.AlignmentFile(args.output, "wb", template=source) as destination:
            for record in source.fetch(until_eof=True):
                match = pattern.search(record.query_name)
                if match is None:
                    raise SystemExit(
                        f"Read name does not match UMI regex: {record.query_name}"
                    )
                record.set_tag("RX", match.group("umi"), value_type="Z", replace=True)
                destination.write(record)
                tagged += 1
    print(f"Tagged {tagged} BAM records with RX")


if __name__ == "__main__":
    main()
