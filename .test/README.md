# SnakeVerse Test Fixture

This directory contains a small RNA-seq-like paired-end fixture used by CI.

The automated tests use `.test/configs/generic.yaml` and
`.test/configs/rnaseq.yaml` as pointer configs. Those configs resolve only files
under `.test/configs/`, so CI can lint and dry-run the workflow without changing
the user-facing `config/` directory.

The workflow currently assumes prebuilt aligner indexes, so CI does not execute
alignment. It validates configs and performs Snakemake dry-runs against the real
FASTQ/FASTA/GTF paths. Full execution tests can be added later once index-building
rules or committed miniature indexes exist.

