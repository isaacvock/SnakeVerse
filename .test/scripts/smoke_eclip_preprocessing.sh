#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 0 ]]; then
    export PATH="$1:$PATH"
fi

python .test/scripts/generate_eclip_fastqs.py
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT

cutadapt \
    -g 'A01=^ATTGCTT' -g 'B06=^ACAAGCC' -e 0 --no-indels \
    -o "$scratch/{name}_R1.fastq.gz" -p "$scratch/{name}_R2.fastq.gz" \
    .test/generated/eclip/encode_pool_R1.fastq.gz \
    .test/generated/eclip/encode_pool_R2.fastq.gz \
    > "$scratch/demux.log"

cat "$scratch/A01_R1.fastq.gz" "$scratch/B06_R1.fastq.gz" > "$scratch/merged_R1.fastq.gz"
cat "$scratch/A01_R2.fastq.gz" "$scratch/B06_R2.fastq.gz" > "$scratch/merged_R2.fastq.gz"

umi_tools extract --extract-method=string --bc-pattern2=NNNNN \
    --stdin="$scratch/merged_R1.fastq.gz" --stdout="$scratch/extracted_R1.fastq.gz" \
    --read2-in="$scratch/merged_R2.fastq.gz" --read2-out="$scratch/extracted_R2.fastq.gz" \
    --log="$scratch/encode_extract.log"

[[ $(gzip -cd "$scratch/extracted_R1.fastq.gz" | sed -n '1p') == *_AAAAA ]]
[[ $(gzip -cd "$scratch/extracted_R2.fastq.gz" | sed -n '2p' | awk '{print length}') -eq 50 ]]

umi_tools extract --extract-method=string --bc-pattern2=NNNNNNNN \
    --stdin=.test/generated/eclip/ecclip_ip_R1.fastq.gz \
    --stdout="$scratch/ecclip_R1.fastq.gz" \
    --read2-in=.test/generated/eclip/ecclip_ip_R2.fastq.gz \
    --read2-out="$scratch/ecclip_R2.fastq.gz" \
    --log="$scratch/ecclip_extract.log"

[[ $(gzip -cd "$scratch/ecclip_R1.fastq.gz" | sed -n '1p') == *_ACGTACGT ]]
[[ $(gzip -cd "$scratch/ecclip_R2.fastq.gz" | sed -n '2p' | awk '{print length}') -eq 0 ]]
