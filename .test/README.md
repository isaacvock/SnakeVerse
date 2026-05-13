# SnakeVerse Test Fixture

This directory contains a small RNA-seq-like paired-end fixture used by CI.

The automated tests use `.test/configs/generic.yaml`,
`.test/configs/rnaseq.yaml`, and `.test/configs/atacseq.yaml` as pointer
configs. Those configs resolve only files under `.test/configs/`, so CI can
lint and dry-run the workflow without changing the user-facing `config/`
directory.

The ATAC-seq fixture under `.test/data/ATACseq/` contains paired-end ENCODE
K562 reads restricted to chr21.

CI validates configs, lints the Snakefile, performs dry-runs, and executes:

- a full single-end generic FASTQ-to-BAM run with Bowtie2 and generated indexes
- a full paired-end RNA-seq run with fastp, STAR, generated indexes,
  transcriptome BAMs, featureCounts standard/strict/full-gene counts, and Salmon
  gene/isoform quantification
- a BWA-MEM2 FASTQ-to-BAM dry-run to exercise the third aligner path
- a full paired-end ATAC-seq run on the K562 ENCODE chr21 fixture with Bowtie2,
  duplicate marking, blacklist/TSS reference fields, MACS3 peaks, and ATAC QC
  target assertions
