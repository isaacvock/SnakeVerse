#!/usr/bin/env bash
set -euo pipefail

source /home/isaac/miniforge3/etc/profile.d/conda.sh
conda activate snakeverse-dev

source_root=$(pwd)
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
mkdir -p "$scratch/repo"
tar --exclude=.git --exclude=.snakemake -cf - -C "$source_root" . \
    | tar -xf - -C "$scratch/repo"
cd "$scratch/repo"
export PYTHONPYCACHEPREFIX=/tmp/snakeverse-pycache

for configfile in "$@"; do
    snakemake --snakefile workflow/Snakefile \
        --configfile "$configfile" \
        --dry-run --cores 2 --printshellcmds
done
