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

python .test/scripts/generate_eclip_fastqs.py
snakemake --snakefile workflow/Snakefile --configfile .test/configs/eclip_encode.yaml \
    --cores 2 --show-failed-logs --latency-wait 30
snakemake --snakefile workflow/Snakefile --configfile .test/configs/ecclip_2026.yaml \
    --cores 2 --show-failed-logs --latency-wait 30
snakemake --snakefile workflow/Snakefile --configfile .test/configs/eclip_preextracted.yaml \
    --cores 2 --show-failed-logs --latency-wait 30

samtools view .test/results/eclip_encode/bam/pre_dedup/PUM2_IP_1.bam | grep -q 'RX:Z:'
test -s .test/results/eclip_encode/tracks/crosslink_sites/PUM2_IP_1.plus.bw
test -s .test/results/eclip_encode/tracks/enrichment/PUM2_IP_1.vs.PUM2_SMI_1.plus.bw
samtools view .test/results/ecclip_2026/bam/pre_dedup/RBP_IP_1.bam | grep -q 'RX:Z:'
test -s .test/results/ecclip_2026/tracks/crosslink_sites/RBP_IP_1.plus.bw
test ! -e .test/results/ecclip_2026/tracks/enrichment/RBP_IP_1.vs.TAG_1.plus.bw
samtools view .test/results/eclip_preextracted/bam/pre_dedup/PRE_IP_1.bam | grep -q 'RX:Z:'
test "$(samtools view -c .test/results/eclip_preextracted/bam/umi_dedup/PRE_IP_1.bam)" \
    -lt "$(samtools view -c .test/results/eclip_preextracted/bam/pre_dedup/PRE_IP_1.bam)"

echo "Full conventional eCLIP, ecCLIP, and pre-extracted integration runs passed."
