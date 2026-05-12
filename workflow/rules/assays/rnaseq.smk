from params import render_tool_params, resource_value, tool_extra
from samples import featurecounts_count_read_pairs, featurecounts_paired_end


RNASEQ_TARGETS = []

if output_enabled("fastq_qc"):
    RNASEQ_TARGETS.extend(fastqc_targets(SAMPLES, RESULTS_DIR))

if output_enabled("trimmed_fastq") and step_enabled("trimming", False):
    RNASEQ_TARGETS.extend(trimmed_fastq_targets(SAMPLES, RESULTS_DIR))

if output_enabled("raw_bam"):
    RNASEQ_TARGETS.extend(raw_bam_targets(SAMPLES, RESULTS_DIR))

if output_enabled("filtered_bam") and step_enabled("bam_filter", True):
    RNASEQ_TARGETS.extend(filtered_bam_targets(SAMPLES, RESULTS_DIR))

if output_enabled("bam_qc") and step_enabled("bam_qc", True):
    RNASEQ_TARGETS.extend(bam_qc_targets(SAMPLES, RESULTS_DIR))

if output_enabled("gene_counts"):
    RNASEQ_TARGETS.append(f"{RESULTS_DIR}/counts/featurecounts/gene_counts.txt")

if output_enabled("transcriptome_bam"):
    RNASEQ_TARGETS.extend(transcriptome_bam_targets(SAMPLES, RESULTS_DIR))

if output_enabled("bigwig") and step_enabled("coverage", False):
    RNASEQ_TARGETS.extend(bigwig_targets(SAMPLES, RESULTS_DIR))

MULTIQC_INPUTS.extend(RNASEQ_TARGETS)

if output_enabled("multiqc"):
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
        str(WORKFLOW_DIR / "envs" / "subread.yaml")
    params:
        annotation=lambda wildcards: config["genome"]["gtf"],
        rendered=lambda wildcards: render_tool_params(
            config,
            "featurecounts",
            overrides={
                "paired_end": featurecounts_paired_end(SAMPLES, config),
                "count_read_pairs": featurecounts_count_read_pairs(SAMPLES, config),
            },
        ),
        extra=lambda wildcards: tool_extra(config, "featurecounts")
    shell:
        """
        mkdir -p $(dirname {output.counts}) $(dirname {log})
        featureCounts -T {threads} -a {params.annotation} -o {output.counts} \
            {params.rendered} {params.extra} {input.bam} > {log} 2>&1
        """
