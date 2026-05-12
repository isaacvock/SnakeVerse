# SnakeVerse Test Fixture

This directory contains a small RNA-seq-like paired-end fixture used by CI.

The automated tests use `.test/configs/generic.yaml` and
`.test/configs/rnaseq.yaml` as pointer configs. Those configs resolve only files
under `.test/configs/`, so CI can lint and dry-run the workflow without changing
the user-facing `config/` directory.

CI validates configs, lints the Snakefile, performs dry-runs, and executes:

- a full single-end generic FASTQ-to-BAM run with Bowtie2 and generated indexes
- a full paired-end RNA-seq run with STAR, generated indexes, transcriptome BAMs,
  and featureCounts gene counts
- a BWA-MEM2 FASTQ-to-BAM dry-run to exercise the third aligner path
