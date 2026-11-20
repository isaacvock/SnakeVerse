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
samtools view .test/results/eclip_encode/bam/pre_dedup/PUM2_IP_1.bam | grep -q 'BC:Z:AAGCAAT'
samtools view .test/results/eclip_encode/bam/pre_dedup/PUM2_IP_1.bam | grep -q 'BC:Z:GGCTTGT'
samtools view .test/results/eclip_encode/bam/pre_dedup/PUM2_IP_1.bam | grep -q '|BC=A01'
test "$(samtools view -c .test/results/eclip_encode/bam/pre_dedup/PUM2_IP_1.bam)" -eq 10
test "$(samtools view -c .test/results/eclip_encode/bam/umi_dedup/PUM2_IP_1.bam)" -eq 8
samtools view .test/results/eclip_encode/bam/umi_dedup/PUM2_IP_1.bam | grep 'RX:Z:AAAAA' | grep -q 'BC:Z:AAGCAAT'
samtools view .test/results/eclip_encode/bam/umi_dedup/PUM2_IP_1.bam | grep 'RX:Z:AAAAA' | grep -q 'BC:Z:GGCTTGT'
test -s .test/results/eclip_encode/tracks/read_coverage/PUM2_IP_1.plus.bw
test -s .test/results/eclip_encode/tracks/crosslink_sites/PUM2_IP_1.plus.bw
test -s .test/results/eclip_encode/tracks/enrichment/PUM2_IP_1.vs.PUM2_SMI_1.plus.bw

samtools view .test/results/ecclip_2026/bam/pre_dedup/RBP_IP_1.bam | grep -q 'RX:Z:'
test "$(samtools view -c .test/results/ecclip_2026/bam/pre_dedup/RBP_IP_1.bam)" -eq 4
test "$(samtools view -c .test/results/ecclip_2026/bam/umi_dedup/RBP_IP_1.bam)" -eq 3
test -s .test/results/ecclip_2026/tracks/read_coverage/RBP_IP_1.plus.bw
test ! -e .test/results/ecclip_2026/tracks/crosslink_sites/RBP_IP_1.plus.bw
test ! -e .test/results/ecclip_2026/tracks/enrichment/RBP_IP_1.vs.TAG_1.plus.bw

samtools view .test/results/eclip_preextracted/bam/pre_dedup/PRE_IP_1.bam | grep -q 'RX:Z:AAAAA'
samtools view .test/results/eclip_preextracted/bam/pre_dedup/PRE_IP_1.bam | grep -q 'BC:Z:AAGCAAT'
samtools view .test/results/eclip_preextracted/bam/pre_dedup/PRE_IP_1.bam | grep -q 'BC:Z:GGCTTGT'
test "$(samtools view -c .test/results/eclip_preextracted/bam/pre_dedup/PRE_IP_1.bam)" -eq 6
test "$(samtools view -c .test/results/eclip_preextracted/bam/umi_dedup/PRE_IP_1.bam)" -eq 4

echo "Full conventional eCLIP, ecCLIP, and pre-extracted integration runs passed."
