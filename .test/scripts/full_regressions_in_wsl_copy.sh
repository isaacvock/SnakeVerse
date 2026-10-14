#!/usr/bin/env bash
set -euo pipefail

source /home/isaac/miniforge3/etc/profile.d/conda.sh
conda activate snakeverse-dev

source_root=$(pwd)
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
mkdir -p "$scratch/repo"
tar --exclude=.git --exclude=.snakemake --exclude=.test/results -cf - -C "$source_root" . \
    | tar -xf - -C "$scratch/repo"
cd "$scratch/repo"
export PYTHONPYCACHEPREFIX=/tmp/snakeverse-pycache

for configfile in \
    .test/configs/generic.yaml \
    .test/configs/rnaseq.yaml \
    .test/configs/atacseq.yaml
do
    snakemake --snakefile workflow/Snakefile --configfile "$configfile" \
        --cores 2 --show-failed-logs --latency-wait 30
done

test -s .test/results/generic/bam/filtered/WT_1.bam
test -s .test/results/rnaseq/counts/featurecounts/gene_counts.txt
test -s .test/results/atacseq/bam/filtered/K562_ATAC_rep1.bam

echo "Generic, RNA-seq, and ATAC-seq regression runs passed."
