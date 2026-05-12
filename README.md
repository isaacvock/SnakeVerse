# SnakeVerse

SnakeVerse is an early-stage, modular Snakemake framework for NGS data processing.
The current scaffold supports two profile stacks:

- generic FASTQ-to-BAM processing for custom or unusual assays
- bulk RNA-seq with STAR alignment and featureCounts gene counting

The repository is intentionally Snakemake-native. SnakeVerse provides workflow
rules, config templates, and a small local config helper, but Snakemake remains
the runner.

## Design Model

One Snakemake run processes one assay/profile stack. Multi-assay projects should
use multiple run configs, for example one RNA-seq run and one ATAC-seq run, rather
than mixing assays in a single sample sheet.

The active entry point is `config/config.yaml`:

```yaml
run_config: config/runs/generic.yaml
```

At runtime, `workflow/lib/config.py` resolves that pointer by loading:

1. the active pointer config
2. the referenced run config
3. each profile in `profile_stack`, in order
4. active tool profiles from `config/profiles/tools/*.yaml`
5. run-level overrides

Later profiles override earlier profiles, and the run config overrides profiles.

## Repository Layout

```text
workflow/
  Snakefile
  rules/
    common/
    assays/
  envs/
  lib/

config/
  config.yaml
  bin/ngsflow.py
  _ngsflow/
    manifest.yaml
    schemas/
    templates/
  runs/
  samples/
  profiles/
```

`workflow/` contains the reusable workflow implementation. `config/_ngsflow/`
contains shipped templates and schemas. `config/runs`, `config/samples`, and
`config/profiles` contain active, user-editable working configs created from
those templates.

This separation matters for Snakedeploy: a deployed local `workflow/` directory
may be only a thin reference to a remote workflow, so the static config helper
and templates live under `config/`.

## Environment

For development in this repository, use the WSL conda environment named
`snakeverse-dev`:

```bash
conda activate snakeverse-dev
```

The workflow itself uses per-rule conda environment YAMLs in `workflow/envs/`.

## Initialize a Generic FASTQ-to-BAM Run

```bash
python config/bin/ngsflow.py list assays
python config/bin/ngsflow.py list presets --assay generic
python config/bin/ngsflow.py init-run \
  --assay generic \
  --preset fastq_to_bam \
  --genome hg38 \
  --run-name generic
```

This copies templates into:

- `config/runs/generic.yaml`
- `config/samples/generic.tsv`
- `config/profiles/assays/`
- `config/profiles/protocols/`
- `config/profiles/tools/`
- `config/profiles/genomes/`

Then edit the sample sheet, genome paths, and tool profiles for your data.

## Initialize an RNA-seq Run

```bash
python config/bin/ngsflow.py init-run \
  --assay rnaseq \
  --preset star_featurecounts \
  --genome hg38 \
  --run-name rnaseq
```

The RNA-seq preset currently includes:

- FastQC
- optional cutadapt trimming
- STAR alignment
- samtools BAM filtering and indexing
- samtools BAM QC
- featureCounts gene-level counts
- optional deepTools BigWig generation
- MultiQC

The RNA-seq sample sheet includes `condition`, `replicate`, and `strandedness`.
The `strandedness` column is available to the workflow for strand-aware coverage
scaffolding. featureCounts strandedness is controlled in the editable
`config/profiles/tools/featurecounts.yaml` profile.

## Activate an Existing Run

```bash
python config/bin/ngsflow.py activate-run rnaseq
```

This updates `config/config.yaml` to point to `config/runs/rnaseq.yaml`.

## Explain and Validate Configs

```bash
python config/bin/ngsflow.py explain --configfile config/config.yaml
python config/bin/ngsflow.py validate --configfile config/config.yaml
```

Validation checks that referenced config files exist, sample sheets have required
columns, relevant tool profiles are active, and genome profiles contain fields
needed by the selected aligner. It also warns about missing FASTQs and reference
paths. Those warnings are expected before you replace template paths with real
local paths.

## Run Snakemake Locally

From the repository root:

```bash
snakemake --configfile config/config.yaml --use-conda --cores 16
```

For a dry-run:

```bash
snakemake --configfile config/config.yaml --dry-run
```

If the template FASTQ or reference paths have not been replaced with real files,
Snakemake may report missing input files. Update `config/samples/<run>.tsv` and
`config/profiles/genomes/<genome>.yaml` before a real run.

## Cluster and Slurm Execution

SnakeVerse does not implement Slurm or cluster submission logic. Use an external
Snakemake workflow profile, such as `smk-simple-slurm`, when you want cluster
execution:

```bash
snakemake \
  --configfile config/config.yaml \
  --use-conda \
  --workflow-profile path/to/smk-simple-slurm
```

Keeping execution profiles external avoids coupling biological workflow logic to
a particular cluster.

## Editing Tool Profiles

Every meaningful tool parameter should be edited in `config/profiles/tools/`,
not buried in a Snakefile. For example:

```yaml
tool: bowtie2
version: "2.5"

params:
  sensitivity: "--very-sensitive"
  max_insert_size: 1000
  no_mixed: false
  no_discordant: false

extra: ""
```

The workflow renders these structured settings into command-line arguments.
Boolean values become flags when true and are omitted when false. Strings and
numbers become flag values. Tool-specific renderers handle common differences
for Bowtie2, STAR, samtools, featureCounts, deepTools, cutadapt, FastQC, and
MultiQC.

`extra` is appended verbatim to the relevant tool command. Use it for flags that
are too new, too specialized, or too awkward to model structurally yet.

## Adding a New Assay or Preset

The intended extension path is explicit rather than magical:

1. Add or edit templates under `config/_ngsflow/templates/`.
2. Register the assay or preset in `config/_ngsflow/manifest.yaml`.
3. Add an assay Snakefile under `workflow/rules/assays/` if new outputs or rules
   are needed.
4. Reuse common rules where possible.
5. Add new per-rule conda envs under `workflow/envs/` only when new tools are
   introduced.
6. Keep execution-profile and cluster behavior outside this repository.

Future assay stacks such as ATAC-seq, eCLIP, Ribo-seq, ChIP-seq/CUT&Tag,
PRO-seq, and SLAM-seq should be expressible as new assay/protocol/tool/genome
profile combinations plus assay-specific rules only where needed.

## Current Limitations

This is first-round infrastructure, not a biologically exhaustive pipeline.

- Sample handling assumes paired-end inputs in the current templates.
- RNA-seq differential expression is not implemented.
- Salmon quantification is not implemented yet.
- Reference preparation and genome index building are not implemented.
- STAR output handling assumes sorted BAM output from the provided STAR profile.
- Validation is intentionally useful but not comprehensive.
- The helper script is optional and local; it is not an installed Python package.

## Quick Smoke-Test Commands

```bash
python config/bin/ngsflow.py list assays
python config/bin/ngsflow.py list presets --assay generic
python config/bin/ngsflow.py init-run --assay generic --preset fastq_to_bam --genome hg38 --run-name generic --overwrite
python config/bin/ngsflow.py explain --configfile config/config.yaml
python config/bin/ngsflow.py validate --configfile config/config.yaml
snakemake --configfile config/config.yaml --dry-run

python config/bin/ngsflow.py init-run --assay rnaseq --preset star_featurecounts --genome hg38 --run-name rnaseq --overwrite
python config/bin/ngsflow.py explain --configfile config/config.yaml
python config/bin/ngsflow.py validate --configfile config/config.yaml
snakemake --configfile config/config.yaml --dry-run
```

