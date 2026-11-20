#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

LABEL_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
SEQ_RE = re.compile(r"^[ACGTRYSWKMBDHVN]+$")
WHITESPACE_SPLIT_RE = re.compile(r"(\s+)")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Stamp inline barcode metadata (|BC=<label>:<sequence>) onto FASTQ "
            "read identifiers and concatenate one or more gzipped FASTQ streams."
        )
    )
    parser.add_argument(
        "--entry",
        action="append",
        required=True,
        help="Input specification formatted as LABEL=SEQUENCE=FASTQ_PATH.",
    )
    parser.add_argument("--output", required=True, help="Output .fastq.gz path.")
    return parser


def parse_entry(raw: str) -> tuple[str, str, Path]:
    parts = raw.split("=", 2)
    if len(parts) != 3 or not all(parts):
        raise SystemExit(
            f"Invalid --entry '{raw}'; expected LABEL=SEQUENCE=FASTQ_PATH"
        )
    label, sequence, path_str = parts[0].strip(), parts[1].strip().upper(), parts[2]
    if not LABEL_RE.fullmatch(label):
        raise SystemExit(f"Invalid barcode label '{label}'")
    if not SEQ_RE.fullmatch(sequence):
        raise SystemExit(f"Invalid barcode sequence '{sequence}' for {label}")
    path = Path(path_str)
    if not path.exists():
        raise SystemExit(f"Input FASTQ does not exist: {path}")
    return label, sequence, path


def stamp_header(line: str, stamp: str) -> str:
    if not line.startswith("@"):
        raise SystemExit(f"Malformed FASTQ header line: {line!r}")
    stripped = line.rstrip("\r\n")
    newline = line[len(stripped) :] or "\n"
    pieces = WHITESPACE_SPLIT_RE.split(stripped[1:], maxsplit=1)
    read_id = pieces[0]
    rest = "".join(pieces[1:])
    return f"@{read_id}{stamp}{rest}{newline}"


def main() -> None:
    args = build_parser().parse_args()
    entries = [parse_entry(item) for item in args.entry]
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    total_records = 0
    with gzip.open(out_path, "wt", encoding="ascii", compresslevel=6) as out_handle:
        for label, sequence, in_path in entries:
            stamp = f"|BC={label}:{sequence}"
            with gzip.open(in_path, "rt", encoding="ascii") as in_handle:
                while True:
                    header = in_handle.readline()
                    if not header:
                        break
                    seq = in_handle.readline()
                    plus = in_handle.readline()
                    qual = in_handle.readline()
                    if not (seq and plus and qual):
                        raise SystemExit(f"Truncated FASTQ record in {in_path}")
                    out_handle.write(stamp_header(header, stamp))
                    out_handle.write(seq)
                    out_handle.write(plus)
                    out_handle.write(qual)
                    total_records += 1
    print(f"Stamped {total_records} FASTQ records into {out_path}")


if __name__ == "__main__":
    main()
