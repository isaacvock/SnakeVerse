from params import render_featurecounts_for_config, render_tool_params, resource_value, tool_extra
from samples import rsem_strandedness_arg, sample_layout


RNASEQ_TARGETS = []

if module_enabled("fastq_qc"):
    RNASEQ_TARGETS.extend(fastqc_targets(SAMPLES, RESULTS_DIR))

if module_enabled("trimming"):
    RNASEQ_TARGETS.extend(trimmed_fastq_targets(SAMPLES, RESULTS_DIR))

if module_enabled("alignment"):
    RNASEQ_TARGETS.extend(raw_bam_targets(SAMPLES, RESULTS_DIR))

if module_enabled("mark_duplicates"):
    RNASEQ_TARGETS.extend(markdup_bam_targets(SAMPLES, RESULTS_DIR))

if module_enabled("bam_filter"):
    RNASEQ_TARGETS.extend(filtered_bam_targets(SAMPLES, RESULTS_DIR))

if module_enabled("bam_qc"):
    RNASEQ_TARGETS.extend(bam_qc_targets(SAMPLES, RESULTS_DIR))

if gene_count_mode_enabled("standard"):
    RNASEQ_TARGETS.append(f"{RESULTS_DIR}/counts/featurecounts/gene_counts.txt")

if gene_count_mode_enabled("exon_strict"):
    RNASEQ_TARGETS.append(f"{RESULTS_DIR}/counts/featurecounts/exon_strict_counts.txt")

if gene_count_mode_enabled("full_gene"):
    RNASEQ_TARGETS.append(f"{RESULTS_DIR}/counts/featurecounts/full_gene_counts.txt")

if module_settings("alignment").get("transcriptome_bam", False):
    RNASEQ_TARGETS.extend(transcriptome_bam_targets(SAMPLES, RESULTS_DIR))

if salmon_result_enabled("isoform_results"):
    RNASEQ_TARGETS.extend(
        f"{RESULTS_DIR}/counts/salmon/{sample}/quant.sf" for sample in SAMPLE_IDS
    )

if salmon_result_enabled("gene_results"):
    RNASEQ_TARGETS.extend(
        f"{RESULTS_DIR}/counts/salmon/{sample}/quant.genes.sf" for sample in SAMPLE_IDS
    )

if rsem_result_enabled("isoform_results"):
    RNASEQ_TARGETS.extend(
        f"{RESULTS_DIR}/counts/rsem/{sample}/{sample}.isoforms.results" for sample in SAMPLE_IDS
    )

if rsem_result_enabled("gene_results"):
    RNASEQ_TARGETS.extend(
        f"{RESULTS_DIR}/counts/rsem/{sample}/{sample}.genes.results" for sample in SAMPLE_IDS
    )

if module_enabled("coverage"):
    RNASEQ_TARGETS.extend(bigwig_targets(SAMPLES, RESULTS_DIR))

MULTIQC_INPUTS.extend(RNASEQ_TARGETS)

if module_enabled("multiqc"):
    RNASEQ_TARGETS.append(f"{RESULTS_DIR}/multiqc/multiqc_report.html")

ASSAY_TARGETS.extend(RNASEQ_TARGETS)


rule featurecounts:
    input:
        bam=lambda wildcards: [final_bam_path(sample) for sample in SAMPLE_IDS],
        bai=lambda wildcards: [final_bai_path(sample) for sample in SAMPLE_IDS]
    output:
        counts=f"{RESULTS_DIR}/counts/featurecounts/gene_counts.txt",
        summary=f"{RESULTS_DIR}/counts/featurecounts/gene_counts.txt.summary"
    log:
        f"{RESULTS_DIR}/logs/featurecounts/gene_counts.log"
    threads:
        int(resource_value(config, "featurecounts", "threads", 8))
    resources:
        mem_mb=int(resource_value(config, "featurecounts", "mem_mb", 8192)),
        runtime_min=int(resource_value(config, "featurecounts", "runtime_min", 180))
    conda:
        workflow_file("envs/subread.yaml")
    params:
        annotation=lambda wildcards: config["reference"]["gtf"],
        rendered=lambda wildcards: render_featurecounts_for_config(config, SAMPLES),
        extra=lambda wildcards: tool_extra(config, "featurecounts")
    shell:
        """
        mkdir -p $(dirname {output.counts}) $(dirname {log})
        featureCounts -T {threads} -a {params.annotation} -o {output.counts} \
            {params.rendered} {params.extra} {input.bam} > {log} 2>&1
        """


rule featurecounts_exon_strict:
    input:
        bam=lambda wildcards: [final_bam_path(sample) for sample in SAMPLE_IDS],
        bai=lambda wildcards: [final_bai_path(sample) for sample in SAMPLE_IDS]
    output:
        counts=f"{RESULTS_DIR}/counts/featurecounts/exon_strict_counts.txt",
        summary=f"{RESULTS_DIR}/counts/featurecounts/exon_strict_counts.txt.summary"
    log:
        f"{RESULTS_DIR}/logs/featurecounts/exon_strict_counts.log"
    threads:
        int(resource_value(config, "featurecounts", "threads", 8))
    resources:
        mem_mb=int(resource_value(config, "featurecounts", "mem_mb", 8192)),
        runtime_min=int(resource_value(config, "featurecounts", "runtime_min", 180))
    conda:
        workflow_file("envs/subread.yaml")
    params:
        annotation=lambda wildcards: config["reference"]["gtf"],
        rendered=lambda wildcards: render_featurecounts_for_config(
            config,
            SAMPLES,
            overrides={"feature_type": "exon", "non_overlap": 0}
        ),
        extra=lambda wildcards: tool_extra(config, "featurecounts")
    shell:
        """
        mkdir -p $(dirname {output.counts}) $(dirname {log})
        featureCounts -T {threads} -a {params.annotation} -o {output.counts} \
            {params.rendered} {params.extra} {input.bam} > {log} 2>&1
        """


rule gene_regions_saf:
    input:
        gtf=lambda wildcards: config["reference"]["gtf"]
    output:
        saf=f"{RESULTS_DIR}/reference/annotation/{GENOME_SLUG}.gene_regions.saf"
    log:
        f"{RESULTS_DIR}/logs/featurecounts/gene_regions_saf.log"
    conda:
        workflow_file("envs/subread.yaml")
    params:
        project_root=PROJECT_ROOT
    script:
        workflow_file("scripts/gtf_to_gene_saf.py")


rule featurecounts_full_gene:
    input:
        annotation=f"{RESULTS_DIR}/reference/annotation/{GENOME_SLUG}.gene_regions.saf",
        bam=lambda wildcards: [final_bam_path(sample) for sample in SAMPLE_IDS],
        bai=lambda wildcards: [final_bai_path(sample) for sample in SAMPLE_IDS]
    output:
        counts=f"{RESULTS_DIR}/counts/featurecounts/full_gene_counts.txt",
        summary=f"{RESULTS_DIR}/counts/featurecounts/full_gene_counts.txt.summary"
    log:
        f"{RESULTS_DIR}/logs/featurecounts/full_gene_counts.log"
    threads:
        int(resource_value(config, "featurecounts", "threads", 8))
    resources:
        mem_mb=int(resource_value(config, "featurecounts", "mem_mb", 8192)),
        runtime_min=int(resource_value(config, "featurecounts", "runtime_min", 180))
    conda:
        workflow_file("envs/subread.yaml")
    params:
        rendered=lambda wildcards: render_featurecounts_for_config(
            config,
            SAMPLES,
            overrides={"annotation_format": "SAF"},
            drop_keys=("feature_type", "attribute_type"),
        ),
        extra=lambda wildcards: tool_extra(config, "featurecounts")
    shell:
        """
        mkdir -p $(dirname {output.counts}) $(dirname {log})
        featureCounts -T {threads} -a {input.annotation} -o {output.counts} \
            {params.rendered} {params.extra} {input.bam} > {log} 2>&1
        """


if salmon_quant_enabled():
    rule salmon_transcripts:
        input:
            fasta=lambda wildcards: config["reference"]["fasta"],
            gtf=lambda wildcards: config["reference"]["gtf"]
        output:
            fasta=f"{RESULTS_DIR}/reference/salmon/{GENOME_SLUG}.transcripts.fa"
        log:
            f"{RESULTS_DIR}/logs/salmon/build_transcripts.log"
        conda:
            workflow_file("envs/salmon.yaml")
        shell:
            """
            mkdir -p $(dirname {output.fasta}) $(dirname {log})
            gffread {input.gtf} -g {input.fasta} -w {output.fasta} > {log} 2>&1
            """

    if salmon_result_enabled("gene_results"):
        rule salmon_quant:
            input:
                bam=f"{RESULTS_DIR}/bam/transcriptome/{{sample}}.bam",
                transcripts=f"{RESULTS_DIR}/reference/salmon/{GENOME_SLUG}.transcripts.fa",
                gtf=lambda wildcards: config["reference"]["gtf"]
            output:
                isoforms=f"{RESULTS_DIR}/counts/salmon/{{sample}}/quant.sf",
                genes=f"{RESULTS_DIR}/counts/salmon/{{sample}}/quant.genes.sf"
            log:
                f"{RESULTS_DIR}/logs/salmon/{{sample}}.log"
            threads:
                int(resource_value(config, "salmon", "threads", 8))
            resources:
                mem_mb=int(resource_value(config, "salmon", "mem_mb", 8192)),
                runtime_min=int(resource_value(config, "salmon", "runtime_min", 180))
            conda:
                workflow_file("envs/salmon.yaml")
            params:
                outdir=lambda wildcards: f"{RESULTS_DIR}/counts/salmon/{wildcards.sample}",
                gene_map=lambda wildcards, input: f"--geneMap {input.gtf}",
                rendered=lambda wildcards: render_tool_params(config, "salmon"),
                extra=lambda wildcards: tool_extra(config, "salmon")
            shell:
                """
                mkdir -p {params.outdir} $(dirname {log})
                salmon quant -p {threads} -t {input.transcripts} -a {input.bam} \
                    {params.rendered} {params.gene_map} {params.extra} \
                    -o {params.outdir} > {log} 2>&1
                """
    else:
        rule salmon_quant:
            input:
                bam=f"{RESULTS_DIR}/bam/transcriptome/{{sample}}.bam",
                transcripts=f"{RESULTS_DIR}/reference/salmon/{GENOME_SLUG}.transcripts.fa"
            output:
                isoforms=f"{RESULTS_DIR}/counts/salmon/{{sample}}/quant.sf"
            log:
                f"{RESULTS_DIR}/logs/salmon/{{sample}}.log"
            threads:
                int(resource_value(config, "salmon", "threads", 8))
            resources:
                mem_mb=int(resource_value(config, "salmon", "mem_mb", 8192)),
                runtime_min=int(resource_value(config, "salmon", "runtime_min", 180))
            conda:
                workflow_file("envs/salmon.yaml")
            params:
                outdir=lambda wildcards: f"{RESULTS_DIR}/counts/salmon/{wildcards.sample}",
                rendered=lambda wildcards: render_tool_params(config, "salmon"),
                extra=lambda wildcards: tool_extra(config, "salmon")
            shell:
                """
                mkdir -p {params.outdir} $(dirname {log})
                salmon quant -p {threads} -t {input.transcripts} -a {input.bam} \
                    {params.rendered} {params.extra} -o {params.outdir} > {log} 2>&1
                """


if rsem_quant_enabled():
    rule rsem_reference:
        input:
            fasta=lambda wildcards: config["reference"]["fasta"],
            gtf=lambda wildcards: config["reference"]["gtf"]
        output:
            marker=touch(f"{RESULTS_DIR}/reference/rsem/{GENOME_SLUG}/.snakeverse_rsem_reference.done")
        log:
            f"{RESULTS_DIR}/logs/rsem/build_reference.log"
        threads:
            int(resource_value(config, "rsem_prepare_reference", "threads", 4))
        resources:
            mem_mb=int(resource_value(config, "rsem_prepare_reference", "mem_mb", 8192)),
            runtime_min=int(resource_value(config, "rsem_prepare_reference", "runtime_min", 180))
        conda:
            workflow_file("envs/rsem.yaml")
        params:
            prefix=lambda wildcards: f"{RESULTS_DIR}/reference/rsem/{GENOME_SLUG}/{GENOME_SLUG}",
            rendered=lambda wildcards: render_tool_params(config, "rsem", section="prepare_reference"),
            extra=lambda wildcards: tool_extra(config, "rsem", section="prepare_reference")
        shell:
            """
            mkdir -p $(dirname {params.prefix}) $(dirname {log})
            rsem-prepare-reference --gtf {input.gtf} {params.rendered} {params.extra} \
                {input.fasta} {params.prefix} > {log} 2>&1
            """

    if rsem_result_enabled("gene_results"):
        rule rsem_quant:
            input:
                bam=f"{RESULTS_DIR}/bam/transcriptome/{{sample}}.bam",
                reference=f"{RESULTS_DIR}/reference/rsem/{GENOME_SLUG}/.snakeverse_rsem_reference.done"
            output:
                genes=f"{RESULTS_DIR}/counts/rsem/{{sample}}/{{sample}}.genes.results",
                isoforms=f"{RESULTS_DIR}/counts/rsem/{{sample}}/{{sample}}.isoforms.results"
            log:
                f"{RESULTS_DIR}/logs/rsem/{{sample}}.log"
            threads:
                int(resource_value(config, "rsem", "threads", 8))
            resources:
                mem_mb=int(resource_value(config, "rsem", "mem_mb", 8192)),
                runtime_min=int(resource_value(config, "rsem", "runtime_min", 180))
            conda:
                workflow_file("envs/rsem.yaml")
            params:
                reference_prefix=lambda wildcards: f"{RESULTS_DIR}/reference/rsem/{GENOME_SLUG}/{GENOME_SLUG}",
                output_prefix=lambda wildcards: f"{RESULTS_DIR}/counts/rsem/{wildcards.sample}/{wildcards.sample}",
                paired=lambda wildcards: "--paired-end" if sample_layout(SAMPLES, wildcards.sample) == "paired" else "",
                strandedness=lambda wildcards: rsem_strandedness_arg(SAMPLES, wildcards.sample),
                rendered=lambda wildcards: render_tool_params(config, "rsem", section="quant"),
                extra=lambda wildcards: tool_extra(config, "rsem", section="quant")
            shell:
                """
                mkdir -p $(dirname {params.output_prefix}) $(dirname {log})
                rsem-calculate-expression --alignments {params.paired} {params.strandedness} -p {threads} \
                    {params.rendered} {params.extra} {input.bam} \
                    {params.reference_prefix} {params.output_prefix} > {log} 2>&1
                """
    else:
        rule rsem_quant:
            input:
                bam=f"{RESULTS_DIR}/bam/transcriptome/{{sample}}.bam",
                reference=f"{RESULTS_DIR}/reference/rsem/{GENOME_SLUG}/.snakeverse_rsem_reference.done"
            output:
                isoforms=f"{RESULTS_DIR}/counts/rsem/{{sample}}/{{sample}}.isoforms.results"
            log:
                f"{RESULTS_DIR}/logs/rsem/{{sample}}.log"
            threads:
                int(resource_value(config, "rsem", "threads", 8))
            resources:
                mem_mb=int(resource_value(config, "rsem", "mem_mb", 8192)),
                runtime_min=int(resource_value(config, "rsem", "runtime_min", 180))
            conda:
                workflow_file("envs/rsem.yaml")
            params:
                reference_prefix=lambda wildcards: f"{RESULTS_DIR}/reference/rsem/{GENOME_SLUG}/{GENOME_SLUG}",
                output_prefix=lambda wildcards: f"{RESULTS_DIR}/counts/rsem/{wildcards.sample}/{wildcards.sample}",
                paired=lambda wildcards: "--paired-end" if sample_layout(SAMPLES, wildcards.sample) == "paired" else "",
                strandedness=lambda wildcards: rsem_strandedness_arg(SAMPLES, wildcards.sample),
                rendered=lambda wildcards: render_tool_params(config, "rsem", section="quant"),
                extra=lambda wildcards: tool_extra(config, "rsem", section="quant")
            shell:
                """
                mkdir -p $(dirname {params.output_prefix}) $(dirname {log})
                rsem-calculate-expression --alignments {params.paired} {params.strandedness} -p {threads} \
                    {params.rendered} {params.extra} {input.bam} \
                    {params.reference_prefix} {params.output_prefix} > {log} 2>&1
                """
