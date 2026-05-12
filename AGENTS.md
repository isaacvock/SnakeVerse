# AGENTS.md

## Project overview

This repository will impelement a modular Snakemake NGS workflow framework. The current planned supported modes are:

- generic FASTQ-to-BAM
- RNA-seq with STAR + featureCounts

The workflow is run with Snakemake. Do not add a custom workflow runner.

## Environment

Use the WSL conda environment named `snakeverse-dev`.

Activate it with:

```bash
conda activate snakeverse-dev