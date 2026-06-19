# SnakeVerse

SnakeVerse is an early-stage, modular Snakemake framework for NGS data processing.
The current scaffold supports three profile stacks:

- generic FASTQ-to-BAM processing for custom or unusual assays
- bulk RNA-seq with STAR alignment and featureCounts gene counting
- ATAC-seq with ENCODE4-oriented Bowtie2 alignment, filtering, MACS3 peaks,
  signal tracks, and QC metrics

The repository is intentionally Snakemake-native. SnakeVerse provides workflow
rules, config templates, and a small local config helper, but Snakemake remains
the runner.

## Quickstart

This example deploys SnakeVerse and initializes a bulk RNA-seq run. The generic
FASTQ-to-BAM and ATAC-seq initialization commands are documented below.

### 1. Install the deployment tools

Install [Miniforge](https://github.com/conda-forge/miniforge) and
[Git](https://git-scm.com/install/) if they are not already available, then
create a small environment containing Snakemake and Snakedeploy:

```bash
mamba create \
  --name deploy_snakemake \
  --channel conda-forge \
  --channel bioconda \
  snakemake snakedeploy

conda activate deploy_snakemake
```

SnakeVerse runs its tools in rule-specific Conda environments. To keep those
environments in one reusable location, rather than separately under every
project, optionally set `SNAKEMAKE_CONDA_PREFIX`:

```bash
export SNAKEMAKE_CONDA_PREFIX="/path/to/snakemake-conda-envs"
mkdir -p "$SNAKEMAKE_CONDA_PREFIX"
```

Add the `export` command to your shell startup file if it should apply to future
sessions. Otherwise, it applies only to the current terminal.

### 2. Deploy SnakeVerse

Create one working directory for this project or dataset, enter it, and deploy
the workflow:

```bash
mkdir -p my-rnaseq-project
cd my-rnaseq-project

conda activate deploy_snakemake
snakedeploy deploy-workflow \
  https://github.com/isaacvock/SnakeVerse.git \
  . \
  --branch main
```

The deployment contains a local `config/` directory for editable project files
and a `workflow/` entry point that references the remote SnakeVerse workflow.

### 3. Initialize a run

Use the local helper to materialize an RNA-seq configuration:

```bash
python config/bin/ngsflow.py init-run \
  --assay rnaseq \
  --preset star_featurecounts \
  --genome hg38 \
  --run-name rnaseq
```

For this example, `ngsflow.py` creates or activates the following files. Edit
them in this priority order:

1. **`config/samples/rnaseq.tsv`**: Replace the example rows with the sample and
   sequencing-unit identifiers, FASTQ paths, condition, replicate, and
   strandedness for the run. For SRA input, provide `sra_id` and `sra_layout`
   and leave the FASTQ columns blank.
2. **`config/profiles/genomes/hg38.yaml`**: Set paths to the reference FASTA,
   GTF, chromosome sizes, and any assay-specific resources. If the selected
   aligner's configured index files do not exist, SnakeVerse builds them at that
   path; a blank index path builds them under the run's results directory.
   RNA-seq does not require the blacklist or TSS BED fields.
3. **`config/runs/rnaseq.yaml`**: Set the project name, results directory, and
   desired outputs. This file also records the assay, preset, sample sheet, and
   ordered profile stack for this run.
4. **`config/profiles/protocols/rnaseq_star_featurecounts.yaml`**: Review enabled
   steps, trimming and alignment tools, strand-aware coverage behavior, and
   per-rule thread, memory, and runtime requests.
5. **`config/profiles/tools/*.yaml`**: Review tool-specific command-line
   behavior. The RNA-seq preset materializes editable profiles for FastQC,
   fastp, cutadapt, SRA Tools, Bowtie2, BWA-MEM2, STAR, samtools,
   featureCounts, Salmon, RSEM, deepTools, and MultiQC. Each tool profile has an
   `extra` field for flags not represented by structured parameters.
6. **`config/profiles/assays/rnaseq.yaml`**: This contains stable assay-level
   defaults and the sample schema. Most runs can leave it unchanged.
7. **`config/config.yaml`**: The helper updates this pointer to
   `config/runs/rnaseq.yaml`. Normally, verify it rather than editing it by
   hand; `activate-run` can switch the pointer later.

The copies under `config/_ngsflow/` are shipped templates. Edit the active files
listed above, not the templates, unless you are developing SnakeVerse itself.

### 4. Validate and run

Inspect the resolved stack, validate obvious path and sample-sheet problems,
then perform a dry-run:

```bash
python config/bin/ngsflow.py explain --configfile config/config.yaml
python config/bin/ngsflow.py validate --configfile config/config.yaml

snakemake \
  --configfile config/config.yaml \
  --use-conda \
  --cores 16 \
  --dry-run
```

Remove `--dry-run` to execute the pipeline:

```bash
snakemake \
  --configfile config/config.yaml \
  --use-conda \
  --cores 16
```

## Quickstart: Sherlock and Slurm HPC

The configuration and initialization steps are the same on Stanford Sherlock.
The differences are where reusable environments are stored, activating the
SnakeVerse base container, and using an external Slurm profile for scalable job
submission.

### 1. Configure persistent environment storage

Your Bash startup file is `$HOME/.bashrc`. Add a shared Conda-environment
location so every compute node can access the same rule environments:

```bash
export SNAKEMAKE_CONDA_PREFIX="$OAK/$USER/snakeverse/conda"
```

Then load the setting and create the directory:

```bash
source "$HOME/.bashrc"
mkdir -p "$SNAKEMAKE_CONDA_PREFIX"
```

You can place the export directly in a Slurm launch script instead if you do
not want it applied to every shell session.

### 2. Deploy and initialize

From a suitable project or group filesystem, repeat the normal deployment and
initialization steps:

```bash
conda activate deploy_snakemake
mkdir -p my-rnaseq-project
cd my-rnaseq-project

snakedeploy deploy-workflow \
  https://github.com/isaacvock/SnakeVerse.git \
  . \
  --branch main

python config/bin/ngsflow.py init-run \
  --assay rnaseq \
  --preset star_featurecounts \
  --genome hg38 \
  --run-name rnaseq
```

Edit the generated samples, genome, run, protocol, and tool profiles in the
same priority order described in the standard quickstart.

### 3. Validate the containerized run

SnakeVerse declares a fixed public base container containing a compatible Linux
user space and Miniforge. Sherlock users should enable both deployment methods:

```bash
snakemake \
  --configfile config/config.yaml \
  --use-apptainer \
  --use-conda \
  --cores 1 \
  --dry-run
```

`--use-apptainer` runs rule jobs in the SnakeVerse base container, while
`--use-conda` creates and activates each rule's tool environment inside that
container. The container image is fixed by the workflow and does not need to be
specified in project configuration.

### 4. Scale out with a Slurm profile

Do not run a full production workflow on a Sherlock login node. Use an external
Snakemake profile and preferably submit the controlling Snakemake process as a
small Slurm job. A profile translates SnakeVerse `threads` and `resources`
values into per-rule `sbatch` submissions.

[`yale_profile`](https://github.com/isaacvock/yale_profile) provides a concrete
collection of Snakemake profile, launch-script, and job-status files for a
Slurm HPC system. Its partition names, modules, email settings, and paths are
Yale-specific, and its resource names must be mapped to SnakeVerse resources
when adapting it for Sherlock. Its current profile also targets Snakemake
versions earlier than 8, so use a compatible Snakemake version or port the same
pattern to a current executor profile. For a maintained alternative, see
[`smk-simple-slurm`](https://github.com/jdblischak/smk-simple-slurm).

With a Sherlock-compatible profile, the launch command has this general form:

```bash
snakemake \
  --configfile config/config.yaml \
  --use-apptainer \
  --use-conda \
  --profile /path/to/sherlock-profile
```

The Slurm profile remains external to SnakeVerse: it controls submission and
cluster policy, while the workflow continues to control biological processing,
software environments, and rule resources.

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

## Core Capabilities

The current basic workflows support:

- paired-end and single-end FASTQ inputs
- optional fastp trimming by default, with cutadapt still selectable
- optional SRA Toolkit FASTQ download from `sra_id` sample-sheet values
- FastQC and MultiQC
- alignment with Bowtie2, STAR, or BWA-MEM2
- automatic Bowtie2, STAR, or BWA-MEM2 index creation when the selected index
  path is left blank
- BAM sorting, indexing, filtering, and samtools QC
- optional deepTools BigWig generation
- RNA-seq gene counts with featureCounts, including optional strict exon and
  full gene-body counting modes
- optional STAR transcriptome-aligned BAMs
- optional Salmon or RSEM gene/isoform quantification from STAR transcriptome
  BAMs, with transcript FASTA/RSEM references generated from the genome FASTA
  and GTF
- ATAC-seq duplicate marking, mitochondrial/blacklist filtering, MACS3
  narrowPeak calls, FRiP, library complexity, fragment length summaries, and an
  optional simple TSS enrichment proxy

SnakeVerse does not perform RNA-seq differential expression analysis. The RNA-seq
workflow stops at aligned BAMs, QC, tracks, and read-count matrices.

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
- optional fastp trimming by default, or cutadapt when `trimming.tool` is set
  to `cutadapt`
- STAR alignment
- samtools BAM filtering and indexing
- samtools BAM QC
- featureCounts gene-level counts
- optional featureCounts strict exon and full gene-body counts
- optional STAR transcriptome-aligned BAMs
- optional Salmon or RSEM gene/isoform quantification from BAM input
- optional deepTools BigWig generation
- MultiQC

The RNA-seq sample sheet includes `condition`, `replicate`, and `strandedness`.
The `strandedness` column is available to the workflow for strand-aware coverage
scaffolding. featureCounts strandedness is controlled in the editable
`config/profiles/tools/featurecounts.yaml` profile.

For single-end data, leave the `fastq_2` value blank in the sample sheet. For
RNA-seq gene counting, all samples in one run should share the same paired-end or
single-end layout because featureCounts is executed once across the run.

Sample sheets may include `sra_id` and `sra_layout`. When `sra_id` is set and
FASTQ paths are blank, SnakeVerse downloads FASTQs with SRA Toolkit; set
`sra_layout: paired` for paired-end SRA accessions.

RNA-seq output switches live in the run config. `gene_counts` is the standard
featureCounts exon count matrix. `exon_strict_counts` adds `--nonOverlap 0`,
and the STAR profile uses `alignEndsType: EndToEnd` by default so soft-clipped
bases do not make otherwise exonic alignments fail that strict overlap test.
`full_gene_counts` counts reads against generated gene-span SAF annotations.
`salmon_gene_quant`, `salmon_isoform_quant`, `rsem_gene_quant`, and
`rsem_isoform_quant` enable transcriptome-BAM quantification.

## Initialize an ATAC-seq Run

```bash
python config/bin/ngsflow.py init-run \
  --assay atacseq \
  --preset encode4_bowtie2_macs3 \
  --genome hg38 \
  --run-name atacseq
```

The ATAC-seq preset follows the practical shape of the ENCODE4 ATAC-seq
processing recommendations: Bowtie2 alignment, high-quality filtered BAMs,
non-mitochondrial/blacklist filtering, normalized signal tracks, peak calls, and
QC metrics such as FRiP, library complexity, fragment length distribution, and
TSS enrichment scaffolding. See the ENCODE ATAC-seq standards and processing
pipeline page for the upstream expectations around biological replicates,
read depth, IDR, FRiP, TSS enrichment, and library complexity:
https://www.encodeproject.org/data-standards/atac-seq/atac-encode4/

The initial SnakeVerse ATAC profile does not yet implement replicated IDR peak
selection or pseudoreplicates. It produces per-sample MACS3 narrowPeak files and
QC outputs that make a future IDR layer straightforward to add.

For single-end ATAC-seq data, leave `fastq_2` blank in the sample sheet. The
default ATAC profile expects two or more biological replicate values in the
`replicate` column and validation will warn when a run has fewer.

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

## Reference Indexes

Genome profiles can point to an existing index or choose where SnakeVerse should
build one. Leave the selected aligner's index path blank to build under the
active results directory:

```yaml
genome:
  name: hg38
  fasta: resources/genomes/hg38/hg38.fa
  gtf: resources/genomes/hg38/gencode.annotation.gtf
  bowtie2_index: ""
  star_index: ""
  bwa_mem2_index: ""
```

If you already have an index, set the corresponding field:

```yaml
genome:
  bowtie2_index: resources/genomes/hg38/bowtie2/hg38
  star_index: resources/genomes/hg38/star
  bwa_mem2_index: resources/genomes/hg38/bwa_mem2/hg38
```

You can also set a path that does not exist yet. SnakeVerse checks for the
aligner's expected index files and, when they are missing, builds the index at
the configured path. For STAR, the configured value is the index directory:

```yaml
genome:
  fasta: resources/genomes/hg38/hg38.fa
  gtf: resources/genomes/hg38/gencode.annotation.gtf
  star_index: resources/genomes/hg38/star
```

STAR index-building uses `genome.fasta` and, when present, `genome.gtf`.
Bowtie2 and BWA-MEM2 index-building use `genome.fasta`.

## Cluster and Slurm Execution

SnakeVerse does not implement Slurm or cluster submission logic. Use an external
Snakemake profile, such as
[`smk-simple-slurm`](https://github.com/jdblischak/smk-simple-slurm), when you
want cluster execution:

```bash
snakemake \
  --configfile config/config.yaml \
  --use-conda \
  --profile path/to/slurm-profile
```

On Sherlock, also pass `--use-apptainer`. The
[`yale_profile`](https://github.com/isaacvock/yale_profile) repository is a
useful example of the profile, controller-job, and status-checking files used to
scale Snakemake across a Slurm cluster; adapt its site-specific settings and
check its stated Snakemake-version compatibility. Keeping execution profiles
external avoids coupling biological workflow logic to a particular cluster.

## Editing Tool Profiles

Every meaningful tool parameter should be edited in `config/profiles/tools/`,
not buried in a Snakefile. For example:

```yaml
tool: bowtie2
version: "2.5"

params:
  index: {}
  align:
    sensitivity: "--very-sensitive"
    max_insert_size: 1000
    no_mixed: false
    no_discordant: false

extra: ""
```

The workflow renders these structured settings into command-line arguments.
Boolean values become flags when true and are omitted when false. Strings and
numbers become flag values. Tool-specific renderers handle common differences
for Bowtie2, STAR, samtools, featureCounts, Salmon, RSEM, deepTools, fastp,
cutadapt, FastQC, MACS3, and MultiQC.

`extra` is appended verbatim to the relevant tool command. Use it for flags that
are too new, too specialized, or too awkward to model structurally yet.

To switch aligners, edit the protocol profile:

```yaml
alignment:
  tool: bwa_mem2
```

Valid values are `bowtie2`, `star`, and `bwa_mem2`. Keep the matching tool
profile in `config/profiles/tools/`.

For STAR transcriptome-aligned BAMs, set `outputs.transcriptome_bam: true`.
Salmon and RSEM quantification outputs also request this BAM internally. In both
cases, make sure the STAR align parameters include `TranscriptomeSAM`:

```yaml
params:
  align:
    quantMode: "GeneCounts TranscriptomeSAM"
```

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

Future assay stacks such as eCLIP, Ribo-seq, ChIP-seq/CUT&Tag, PRO-seq, and
SLAM-seq should be expressible as new assay/protocol/tool/genome profile
combinations plus assay-specific rules only where needed.

## Current Limitations

This is first-round infrastructure, not a biologically exhaustive pipeline.

- RNA-seq differential expression is not implemented.
- ATAC-seq replicated IDR, pseudoreplicate peak selection, and exact ENCODE
  TSS enrichment scoring are not implemented yet.
- Reference FASTA/GTF acquisition is not implemented.
- STAR output handling assumes sorted genome BAM output from the provided STAR profile.
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

python config/bin/ngsflow.py init-run --assay atacseq --preset encode4_bowtie2_macs3 --genome hg38 --run-name atacseq --overwrite
python config/bin/ngsflow.py explain --configfile config/config.yaml
python config/bin/ngsflow.py validate --configfile config/config.yaml
snakemake --configfile config/config.yaml --dry-run
```
