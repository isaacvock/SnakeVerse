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

echo "=== 1. Python syntax & pytest ==="
python -m py_compile \
    config/bin/ngsflow.py \
    workflow/lib/config.py \
    workflow/lib/params.py \
    workflow/lib/refs.py \
    workflow/lib/samples.py \
    workflow/lib/clip.py \
    workflow/lib/mutprof.py \
    workflow/lib/tracks.py \
    workflow/lib/validation.py \
    workflow/scripts/call_snps.py \
    workflow/scripts/count_mutations.py \
    workflow/scripts/build_editing_tracks.py
pytest -q tests/test_eclip.py tests/test_mutational_profiling.py

echo "=== 2. Generate synthetic fixtures & validate configs ==="
python .test/scripts/generate_eclip_fastqs.py
python .test/scripts/generate_mutprof_fastqs.py

python config/bin/ngsflow.py list assays
python config/bin/ngsflow.py list recipes --assay mutational_profiling
python config/bin/ngsflow.py explain mutations --configfile .test/configs/mutprof_dual_controls.yaml
python config/bin/ngsflow.py validate --configfile .test/configs/rnaseq.yaml
python config/bin/ngsflow.py validate --configfile .test/configs/eclip_encode.yaml
python config/bin/ngsflow.py validate --configfile .test/configs/mutprof_dual_controls.yaml
python config/bin/ngsflow.py validate --configfile .test/configs/mutprof_unstranded_nocontrols.yaml

echo "=== 3. Recipe initialization smoke test ==="
smoke_dir=$(mktemp -d)
cp -a config "$smoke_dir/config"
(
    cd "$smoke_dir"
    python config/bin/ngsflow.py init-run --assay mutational_profiling --recipe hypertribe --genome hg38 --run-name smoke-ht
    python config/bin/ngsflow.py init-run --assay mutational_profiling --recipe stamp --genome hg38 --run-name smoke-stamp
    python config/bin/ngsflow.py init-run --assay mutational_profiling --recipe dual_deaminase --genome hg38 --run-name smoke-dual
    python config/bin/ngsflow.py explain
    python config/bin/ngsflow.py explain mutations
)
rm -rf "$smoke_dir"

echo "=== 4. Snakemake lint ==="
snakemake --snakefile workflow/Snakefile --configfile .test/configs/rnaseq.yaml --lint
snakemake --snakefile workflow/Snakefile --configfile .test/configs/mutprof_dual_controls.yaml --lint
snakemake --snakefile workflow/Snakefile --configfile .test/configs/mutprof_unstranded_nocontrols.yaml --lint

echo "=== 5. Snakemake dry-runs ==="
snakemake --snakefile workflow/Snakefile --configfile .test/configs/rnaseq.yaml --dry-run --cores 2
snakemake --snakefile workflow/Snakefile --configfile .test/configs/mutprof_dual_controls.yaml --dry-run --cores 2
snakemake --snakefile workflow/Snakefile --configfile .test/configs/mutprof_unstranded_nocontrols.yaml --dry-run --cores 2

echo "=== 6. Full end-to-end mutational_profiling integration runs ==="
snakemake --snakefile workflow/Snakefile --configfile .test/configs/mutprof_dual_controls.yaml \
    --cores 2 --show-failed-logs --latency-wait 30
snakemake --snakefile workflow/Snakefile --configfile .test/configs/mutprof_unstranded_nocontrols.yaml \
    --cores 2 --show-failed-logs --latency-wait 30

python - <<'PY'
import csv
import gzip
from pathlib import Path
import pyBigWig
import pysam

with pysam.AlignmentFile(".test/results/mutprof_dual_controls/bam/raw/RBP_1.bam", "rb") as sf:
    assert sum(1 for _ in sf) == 32
with pysam.AlignmentFile(".test/results/mutprof_dual_controls/bam/filtered/RBP_1.bam", "rb") as sf:
    assert sum(1 for _ in sf) == 28

snps_ctrl = Path(".test/results/mutprof_dual_controls/mutations/anomalies/snps.txt").read_text().splitlines()
assert "A:G:chrMut:1156" in snps_ctrl
assert "A:G:chrMut:1151" not in snps_ctrl

with gzip.open(".test/results/mutprof_dual_controls/mutations/per_base/RBP_1.tsv.gz", "rt") as fh:
    rows = {(int(r["pos"]), r["strand"]): r for r in csv.DictReader(fh, delimiter="\t")}

assert (1156, "+") not in rows, "Germline SNP pos 1156 must be masked"
p1151 = rows[(1151, "+")]
assert p1151["ref_base"] == "A" and int(p1151["depth"]) == 8 and int(p1151["A"]) == 4 and int(p1151["G"]) == 4

p5151 = rows[(5151, "-")]
assert p5151["ref_base"] == "C" and int(p5151["depth"]) == 6 and int(p5151["C"]) == 3 and int(p5151["T"]) == 3

with pyBigWig.open(".test/results/mutprof_dual_controls/tracks/editing/RBP_1.AG.edit_rate.plus.bw") as bw:
    assert abs(bw.values("chrMut", 1150, 1151)[0] - 0.5) < 1e-5
with pyBigWig.open(".test/results/mutprof_dual_controls/tracks/editing/RBP_1.CT.edit_count.minus.bw") as bw:
    assert abs(bw.values("chrMut", 5150, 5151)[0] - 3.0) < 1e-5
assert Path(".test/results/mutprof_dual_controls/tracks/editing/RBP_1.AG.edit_rate.plus.bedGraph.gz").exists()
assert Path(".test/results/mutprof_dual_controls/counts/featurecounts/gene_counts.txt").exists()

snps_unstr = set(
    Path(".test/results/mutprof_unstranded_nocontrols/mutations/anomalies/snps.txt").read_text().splitlines()
)
assert "A:G:chrMut:10136" in snps_unstr, "100% homozygous AG SNP must be masked"
assert "C:A:chrMut:10141" in snps_unstr, "50% non-target C->A SNP must be masked"
assert "A:G:chrMut:10131" not in snps_unstr, "50% true AG editing site must NOT be masked"

with gzip.open(".test/results/mutprof_unstranded_nocontrols/mutations/per_base/UNSTR_RBP_1.tsv.gz", "rt") as fh:
    u_rows = {(int(r["pos"]), r["strand"]): r for r in csv.DictReader(fh, delimiter="\t")}

assert (10136, ".") not in u_rows
assert (10141, ".") not in u_rows
u10131 = u_rows[(10131, ".")]
assert u10131["ref_base"] == "A" and int(u10131["depth"]) == 10 and int(u10131["A"]) == 5 and int(u10131["G"]) == 5

with pyBigWig.open(".test/results/mutprof_unstranded_nocontrols/tracks/editing/UNSTR_RBP_1.AG.edit_rate.plus.bw") as bw:
    assert abs(bw.values("chrMut", 10130, 10131)[0] - 0.5) < 1e-5

print("All mutational_profiling end-to-end integration assertions passed!")
PY
