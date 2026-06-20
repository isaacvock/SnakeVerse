# SnakeVerse

SnakeVerse is a modular, Snakemake-native framework for common NGS processing.
It currently supports:

- generic FASTQ-to-BAM processing with Bowtie2, STAR, or BWA-MEM2
- bulk RNA-seq with STAR and featureCounts, plus optional Salmon or RSEM
- ENCODE-oriented ATAC-seq with Bowtie2, MACS3, signal tracks, and QC

Snakemake remains the runner. The local `ngsflow.py` script only creates,
explains, and validates configuration files; it never submits or executes the
workflow.

## Configuration in One Minute

One Snakemake invocation processes one assay. A project containing RNA-seq and
ATAC-seq should have two run files and two Snakemake invocations, not a mixed
sample sheet.

SnakeVerse configuration has four active, editable parts:

```text
config/config.yaml          active-run pointer
config/runs/<run>.yaml      complete pipeline and resource configuration
config/samples/<run>.tsv    samples and input reads
config/references/<ref>.yaml reference files and indexes
config/tools/<tool>.yaml    available advanced tool settings
```

`config/_catalog/` is the shipped catalog of recipes, schemas, reference
templates, and tool defaults. Its leading underscore is intentional: ordinary
users should not need to edit it. `init-run` copies the recipe plus an
assay-compatible toolbox into the active directories above.

### Recipes Are Generators

A recipe creates a complete run file. After initialization, the recipe has no
runtime authority and there is no hidden inheritance or merge order. The
`created_from` block records provenance only.

This is a breaking replacement for the earlier `profile_stack`, preset,
`steps`, and `outputs` model. For an existing deployment, initialize a fresh
run and transfer sample paths, reference paths, and intentional tool settings;
old layered run files are rejected rather than ambiguously reinterpreted.

The run's `modules` block is the sole place that enables pipeline behavior:

```yaml
modules:
  trimming:
    enabled: true
    tool: fastp
    adapters:
      mode: auto
      read1:
      read2:

  alignment:
    enabled: true
    tool: star

  gene_counts:
    enabled: true
    tool: featurecounts
```

There are no separate `steps` and `outputs` switches to keep synchronized.
Enabled modules produce their canonical outputs, and Snakemake infers the
upstream work required to make them.

## Quickstart

### 1. Install Snakemake

Install [Miniforge](https://github.com/conda-forge/miniforge), then create a
small deployment environment:

```bash
mamba create -n deploy_snakemake -c conda-forge -c bioconda \
  snakemake snakedeploy pyyaml
conda activate deploy_snakemake
```

For a direct clone, Snakedeploy is optional. For a deployed project:

```bash
mkdir my-rnaseq-project
cd my-rnaseq-project
snakedeploy deploy-workflow \
  https://github.com/isaacvock/SnakeVerse.git \
  . \
  --branch main
```

### 2. Initialize a Run

List the available choices and create an RNA-seq run:

```bash
python config/bin/ngsflow.py list assays
python config/bin/ngsflow.py list recipes --assay rnaseq

python config/bin/ngsflow.py init-run \
  --assay rnaseq \
  --recipe star_featurecounts \
  --genome hg38 \
  --run-name rnaseq
```

Initialization creates or activates these files. Edit them in this order:

1. **`config/samples/rnaseq.tsv`**: Replace every example row. Set sample and
   unit IDs, FASTQ paths, condition, replicate, and strandedness.
2. **`config/references/hg38.yaml`**: Set the FASTA and annotation paths needed
   by this run. Existing index paths are optional; SnakeVerse can build them.
3. **`config/runs/rnaseq.yaml`**: Review adapters first, then enabled modules,
   project/result names, quantification behavior, and resource requests.
4. **`config/tools/star.yaml`**: Review alignment policy, especially read
   length-dependent index settings, multimapping, mismatch, and end alignment.
5. **Other files in `config/tools/`**: Change advanced tool behavior only when
   its defaults are inappropriate for your data.
6. **`config/config.yaml`**: Normally do not edit this by hand. It simply points
   to `config/runs/rnaseq.yaml` and can be changed with `activate-run`.

Initialization copies the assay's compatible toolbox, including optional
trimmers, aligners, downloaders, and quantifiers. A tool file being present
means the tool is available; only the run's modules and sample inputs determine
which tools are actually used.

Existing shared reference and tool files are preserved automatically when
another run is initialized. Existing run and sample files still require
`--skip-existing` to keep them or `--overwrite` to replace them. Because
`--overwrite` also restores shared files to catalog defaults, use it carefully.

### 3. Understand and Validate It

The default explanation is an editing-oriented overview:

```bash
python config/bin/ngsflow.py explain
python config/bin/ngsflow.py validate
```

Ask focused questions when needed:

```bash
python config/bin/ngsflow.py explain trimming
python config/bin/ngsflow.py explain alignment
python config/bin/ngsflow.py explain reference
python config/bin/ngsflow.py explain tool star
```

`explain` reports which file owns each setting and separates tools used by the
run from available but inactive tools. `validate` checks the run, sample
columns, selected tools, adapter policy, contextual reference requirements,
layouts, strandedness, and obvious missing paths.

`explain tool <name>` prints the active tool file with its comments intact.
Those comments describe each shipped parameter, the arguments SnakeVerse
already supplies, and a link to the tool's own documentation. For example:

```bash
python config/bin/ngsflow.py explain tool rsem
python config/bin/ngsflow.py explain tool star
```

### 4. Run Snakemake

Dry-run first, then remove `--dry-run`:

```bash
snakemake \
  --configfile config/config.yaml \
  --use-conda \
  --cores 16 \
  --dry-run

snakemake \
  --configfile config/config.yaml \
  --use-conda \
  --cores 16
```

## Other Recipes

Generic FASTQ-to-BAM:

```bash
python config/bin/ngsflow.py init-run \
  --assay generic \
  --recipe fastq_to_bam \
  --genome hg38 \
  --run-name generic
```

ATAC-seq:

```bash
python config/bin/ngsflow.py init-run \
  --assay atacseq \
  --recipe encode4_bowtie2_macs3 \
  --genome hg38 \
  --run-name atacseq
```

Switch among initialized runs with:

```bash
python config/bin/ngsflow.py activate-run rnaseq
```

## Setting Ownership

The run file owns settings that define biological intent or pipeline shape:

- enabled modules and selected tools
- adapter mode and adapter sequences
- BAM inclusion, exclusion, duplicate, and contig policy
- RNA counting feature, attribute, strandedness, and paired-read semantics
- transcriptome BAM and transcript quantification outputs
- per-rule threads, memory, and runtime requests

Tool files own advanced command behavior:

- STAR splice alignment, multimapping, mismatch, and output parameters
- Bowtie2 sensitivity and concordance policy
- fastp or cutadapt quality and length filters
- deepTools normalization
- MACS3 peak-calling details
- other tool-specific switches not central to run intent

This boundary is enforced for high-risk duplicate settings. For example,
adapter sequences in `fastp.yaml` or strandedness in `featurecounts.yaml` are
validation errors because their authoritative locations are in the run file.

Tool files are intended to be readable without consulting workflow code. Their
comments identify workflow-controlled arguments such as paths, threads,
pairedness, and strandedness. Under `params`, `false`, `null`, and blank values
are omitted; `true` emits a flag; and strings or numbers emit an option with a
value. Tool-specific renderers translate friendly names such as
`minimum_overlap` to the CLI's actual spelling.

Every tool file also includes `extra` for specialized flags not represented by
structured parameters. A single-command tool uses a string. A multi-command
tool uses named strings so an option reaches only the intended command:

```yaml
tool: bowtie2
params:
  align:
    sensitivity: --very-sensitive
extra:
  index: ""
  align: "--no-unal"
```

Prefer documented `params` entries when available. `extra` is deliberately an
escape hatch: SnakeVerse appends it verbatim and cannot validate the flags.
RSEM, STAR, Bowtie2, BWA-MEM2, samtools, and SRA Toolkit profiles label their
command stages explicitly where needed.

## Adapter Trimming

Adapter policy is deliberately prominent:

```yaml
modules:
  trimming:
    enabled: true
    tool: fastp
    adapters:
      mode: auto       # auto, sequences, or none
      read1:
      read2:
```

`auto` uses fastp detection. `sequences` requires `read1` and also `read2` for
paired-end input. `none` retains quality/length filtering but skips adapter
removal. Automatic detection is not available with cutadapt, so provide
sequences when selecting that tool. The active fastp and cutadapt files contain
quality and minimum-length defaults that should also be reviewed.

## Samples and Read Layout

Single-end and paired-end inputs are supported. Leave `fastq_2` blank for
single-end data. Multiple units for one sample may be supplied as multiple rows,
but units belonging to one sample must share a layout.

RNA-seq adds `condition`, `replicate`, and `strandedness`. Valid strandedness
values are `unstranded`, `forward`, and `reverse`. With
`modules.gene_counts.strandedness: sample`, the sample sheet controls the
featureCounts `-s` value. Because featureCounts creates one matrix per run, all
samples in that run must share a layout and strandedness.

SRA accessions are detected per sample row. Set `sra_id`, leave the FASTQ paths
blank, and set `sra_layout` to `single` or `paired`. The assay toolbox already
contains the SRA tool configuration.

Rows with local FASTQ paths continue to use those files, so local and SRA-backed
samples can coexist in one run. Downloads use fasterq-dump plus parallel pigz
compression.

## References and Indexes

Reference templates are intentionally honest: unavailable files start blank.
Requirements depend on the enabled modules.

- FASTA is required when the selected alignment index must be built.
- GTF is required for RNA featureCounts and transcript quantification.
- chromosome sizes and TSS BED are required for ATAC TSS enrichment.
- a blacklist is optional generally and recommended for ATAC-seq.

```yaml
reference:
  name: hg38
  fasta: resources/hg38/hg38.fa
  gtf: resources/hg38/genes.gtf
  chrom_sizes:
  blacklist:
  tss_bed:
  indexes:
    bowtie2:
    star: resources/hg38/star
    bwa_mem2:
  effective_genome_size: 2913022398
  macs3_genome_size: hs
```

An index path may be blank, ready, or a desired destination that does not yet
exist. SnakeVerse checks for expected index files. If they are absent, it builds
the index at the configured path; a blank path builds under the run's results
directory. STAR uses the FASTA and optional GTF. Bowtie2 and BWA-MEM2 use the
FASTA.

## RNA-seq Behavior

The default RNA-seq recipe produces FastQC, trimmed reads, STAR genome BAMs,
filtered/indexed BAMs, samtools QC, featureCounts gene counts, BigWigs, and
MultiQC. It can also produce:

- strict exon counts
- counts across generated full gene spans
- STAR transcript-coordinate BAMs
- Salmon or RSEM gene and isoform quantification

Enable these in the matching run modules; their tool files are already
available. STAR `quantMode` must include `TranscriptomeSAM` for transcriptome
BAM, Salmon, or RSEM output. SnakeVerse intentionally stops at read counts and
does not perform differential expression analysis.

## ATAC-seq Behavior

The ATAC recipe follows the practical shape of ENCODE4 guidance: paired or
single-end alignment, duplicate marking, mitochondrial and blacklist filtering,
signal tracks, per-sample MACS3 narrowPeak calls, FRiP, library complexity,
fragment lengths, and optional TSS enrichment. See the [ENCODE ATAC-seq
standards](https://www.encodeproject.org/data-standards/atac-seq/atac-encode4/).

Replicated IDR, pseudoreplicates, and exact ENCODE TSS scoring are not yet
implemented.

## Sherlock and Apptainer

Sherlock has an older host glibc, so SnakeVerse declares a fixed base container
with a modern Linux user space and Miniforge. General users can use Conda alone;
Sherlock users should enable both Apptainer and Conda.

Your Bash startup file is `$HOME/.bashrc`. Add reusable Conda and Apptainer
locations on shared/scratch storage:

```bash
export SNAKEMAKE_CONDA_PREFIX="$OAK/$USER/snakeverse/conda"
export APPTAINER_CACHEDIR="$SCRATCH/.apptainer/cache"
```

Load the settings and create both directories:

```bash
source "$HOME/.bashrc"
mkdir -p "$SNAKEMAKE_CONDA_PREFIX" "$APPTAINER_CACHEDIR"
```

Then validate the deployment with:

```bash
snakemake \
  --configfile config/config.yaml \
  --use-apptainer \
  --use-conda \
  --cores 1 \
  --dry-run
```

`--use-apptainer` supplies the modern base user space. `--use-conda` creates and
activates each rule's exact tool environment inside it. The image is fixed in
the workflow and is not a user configuration option.

For production, submit through an external Snakemake Slurm profile. SnakeVerse
does not contain cluster submission policy. The
[`yale_profile`](https://github.com/isaacvock/yale_profile) repository shows the
controller, profile, and status-check structure used to deploy Snakemake jobs
on a Slurm HPC system such as Sherlock. Adapt its site-specific partitions,
paths, modules, resource mapping, and documented Snakemake version. A maintained
alternative is [`smk-simple-slurm`](https://github.com/jdblischak/smk-simple-slurm).

```bash
snakemake \
  --configfile config/config.yaml \
  --use-apptainer \
  --use-conda \
  --profile /path/to/sherlock-profile
```

## Repository Layout

```text
workflow/
  Snakefile
  rules/common/
  rules/assays/
  envs/
  lib/

config/
  config.yaml
  bin/ngsflow.py
  _catalog/
    manifest.yaml
    recipes/
    samples/
    references/
    tools/
    schemas/
  runs/
  samples/
  references/
  tools/
```

Keeping the static helper and catalog under `config/` is important for
Snakedeploy, where the local `workflow/` may only be a thin reference to a
remote source.

## Adding an Assay or Recipe

1. Add a complete run template to `config/_catalog/recipes/`.
2. Add its sample template and any new tool or reference templates.
3. Register the recipe and assay-compatible toolbox in `config/_catalog/manifest.yaml`.
4. Reuse common rules, adding `workflow/rules/assays/<assay>.smk` only for new
   biological outputs.
5. Add a rule-level Conda environment when introducing a new tool.
6. Add validation for new contextual requirements and CI fixtures for the path.

Execution profiles remain external. Future eCLIP, Ribo-seq, ChIP-seq/CUT&Tag,
PRO-seq, and SLAM-seq support should follow the same recipe-to-complete-run
pattern.

## Current Limitations

- Reference FASTA, annotation, blacklist, and TSS resource acquisition is not
  automated.
- RNA-seq differential expression is intentionally out of scope.
- ATAC replicated IDR and pseudoreplicate processing are not implemented.
- Validation catches common configuration errors but is not a substitute for
  reviewing assay and tool parameters.
- `ngsflow.py` is a local helper, not an installed package or workflow runner.
