# SnakeVerse Test Fixture

This directory contains a small RNA-seq-like paired-end fixture used by CI.

The automated tests use `.test/configs/generic.yaml`,
`.test/configs/rnaseq.yaml`, and `.test/configs/atacseq.yaml` as pointer
configs. Those configs resolve only files under `.test/configs/`, so CI can
lint and dry-run the workflow without changing the user-facing `config/`
directory.

CI validates configs, lints the Snakefile, performs dry-runs, and executes:

- a full single-end generic FASTQ-to-BAM run with Bowtie2 and generated indexes
- a full paired-end RNA-seq run with STAR, generated indexes, transcriptome BAMs,
  and featureCounts gene counts
- a BWA-MEM2 FASTQ-to-BAM dry-run to exercise the third aligner path
- an ATAC-seq validate/lint/dry-run path with Bowtie2, duplicate marking,
  blacklist/TSS reference fields, MACS2 peaks, and ATAC QC targets
