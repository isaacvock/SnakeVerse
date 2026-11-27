from mutprof import (
    anomaly_bam_samples,
    editing_track_patterns,
    render_anomaly_calling_args,
    render_editing_tracks_args,
    render_mutation_counts_args,
)
from params import resource_value
from samples import sample_strandedness


if module_enabled("anomaly_calling") or module_enabled("editing_tracks"):
    rule index_reference_fasta:
        input:
            fasta=lambda wildcards: config["reference"]["fasta"]
        output:
            fasta=f"{RESULTS_DIR}/reference/fasta/{GENOME_SLUG}.fa",
            fai=f"{RESULTS_DIR}/reference/fasta/{GENOME_SLUG}.fa.fai"
        log:
            f"{RESULTS_DIR}/logs/samtools/faidx.{GENOME_SLUG}.log"
        threads: 1
        resources:
            mem_mb=int(resource_value(config, "samtools_faidx", "mem_mb", 2048)),
            runtime_min=int(resource_value(config, "samtools_faidx", "runtime_min", 30))
        conda:
            workflow_file("envs/samtools.yaml")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.fasta:q}) $(dirname {log:q})
            ln -sf "$(realpath {input.fasta:q})" {output.fasta:q}
            samtools faidx {output.fasta:q} > {log:q} 2>&1
            """


if module_enabled("anomaly_calling"):
    rule call_anomalies:
        input:
            bams=lambda wildcards: [
                final_bam_path(sample) for sample in anomaly_bam_samples(config, SAMPLES)
            ],
            bais=lambda wildcards: [
                final_bai_path(sample) for sample in anomaly_bam_samples(config, SAMPLES)
            ],
            fasta=f"{RESULTS_DIR}/reference/fasta/{GENOME_SLUG}.fa",
            fai=f"{RESULTS_DIR}/reference/fasta/{GENOME_SLUG}.fa.fai"
        output:
            snps=f"{RESULTS_DIR}/mutations/anomalies/snps.txt",
            tsv=f"{RESULTS_DIR}/mutations/anomalies/anomalies.tsv"
        log:
            f"{RESULTS_DIR}/logs/mutational_profiling/call_anomalies.log"
        threads:
            int(resource_value(config, "anomaly_calling", "threads", 4))
        resources:
            mem_mb=int(resource_value(config, "anomaly_calling", "mem_mb", 8192)),
            runtime_min=int(resource_value(config, "anomaly_calling", "runtime_min", 180))
        conda:
            workflow_file("envs/mutational_profiling.yaml")
        params:
            rendered=lambda wildcards: render_anomaly_calling_args(config, SAMPLES),
            script=lambda wildcards: workflow_file("scripts/call_snps.py")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.snps:q}) $(dirname {log:q})
            python {params.script:q} \
                -f {input.fasta:q} \
                --fai {input.fai:q} \
                -o {output.snps:q} \
                --output-tsv {output.tsv:q} \
                -t {threads} \
                {params.rendered} \
                -b {input.bams:q} > {log:q} 2>&1
            """


if module_enabled("mutation_counts"):
    rule sort_bam_by_qname_for_mutations:
        input:
            bam=lambda wildcards: final_bam_path(wildcards.sample),
            bai=lambda wildcards: final_bai_path(wildcards.sample)
        output:
            bam=temp(f"{RESULTS_DIR}/mutations/tmp/{{sample}}.qname.bam")
        log:
            f"{RESULTS_DIR}/logs/mutational_profiling/qname_sort.{{sample}}.log"
        threads:
            int(resource_value(config, "samtools_qname_sort", "threads", 4))
        resources:
            mem_mb=int(resource_value(config, "samtools_qname_sort", "mem_mb", 8192)),
            runtime_min=int(resource_value(config, "samtools_qname_sort", "runtime_min", 120))
        conda:
            workflow_file("envs/samtools.yaml")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
            samtools sort -@ {threads} -n -o {output.bam:q} {input.bam:q} > {log:q} 2>&1
            """

    rule count_mutations:
        input:
            bam=f"{RESULTS_DIR}/mutations/tmp/{{sample}}.qname.bam",
            snps=(
                [f"{RESULTS_DIR}/mutations/anomalies/snps.txt"]
                if module_enabled("anomaly_calling")
                else []
            )
        output:
            perbase=f"{RESULTS_DIR}/mutations/per_base/{{sample}}.tsv.gz"
        log:
            f"{RESULTS_DIR}/logs/mutational_profiling/count_mutations.{{sample}}.log"
        threads:
            int(resource_value(config, "mutation_counts", "threads", 4))
        resources:
            mem_mb=int(resource_value(config, "mutation_counts", "mem_mb", 8192)),
            runtime_min=int(resource_value(config, "mutation_counts", "runtime_min", 240))
        conda:
            workflow_file("envs/mutational_profiling.yaml")
        params:
            strandedness=lambda wildcards: sample_strandedness(SAMPLES, wildcards.sample),
            snp_arg=lambda wildcards, input: (
                f"--snps {input.snps[0]}" if len(input.snps) > 0 else ""
            ),
            rendered=lambda wildcards: render_mutation_counts_args(config),
            script=lambda wildcards: workflow_file("scripts/count_mutations.py")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.perbase:q}) $(dirname {log:q})
            python {params.script:q} \
                -b {input.bam:q} \
                -o {output.perbase:q} \
                -t {threads} \
                --strandedness {params.strandedness:q} \
                {params.snp_arg} \
                {params.rendered} > {log:q} 2>&1
            """


if module_enabled("editing_tracks"):
    rule build_editing_tracks:
        input:
            perbase=f"{RESULTS_DIR}/mutations/per_base/{{sample}}.tsv.gz",
            fai=f"{RESULTS_DIR}/reference/fasta/{GENOME_SLUG}.fa.fai"
        output:
            tracks=editing_track_patterns(config, RESULTS_DIR)
        wildcard_constraints:
            sample="[^/]+"
        log:
            f"{RESULTS_DIR}/logs/mutational_profiling/editing_tracks.{{sample}}.log"
        threads: 1
        resources:
            mem_mb=int(resource_value(config, "editing_tracks", "mem_mb", 4096)),
            runtime_min=int(resource_value(config, "editing_tracks", "runtime_min", 60))
        conda:
            workflow_file("envs/mutational_profiling.yaml")
        params:
            strandedness=lambda wildcards: sample_strandedness(SAMPLES, wildcards.sample),
            rendered=lambda wildcards: render_editing_tracks_args(config),
            script=lambda wildcards: workflow_file("scripts/build_editing_tracks.py")
        shell:
            """
            set -euo pipefail
            outdir=$(dirname {output.tracks[0]:q})
            mkdir -p "$outdir" $(dirname {log:q})
            python {params.script:q} \
                --perbase {input.perbase:q} \
                --fai {input.fai:q} \
                --strandedness {params.strandedness:q} \
                --out-dir "$outdir" \
                --sample {wildcards.sample:q} \
                {params.rendered} > {log:q} 2>&1
            """
