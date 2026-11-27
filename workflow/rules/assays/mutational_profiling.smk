from mutprof import anomaly_targets, editing_track_targets, per_base_targets
from samples import (
    bam_qc_targets,
    bigwig_targets,
    filtered_bam_targets,
    input_fastqc_targets,
    markdup_bam_targets,
    pre_dedup_bam_targets,
    pre_dedup_qc_targets,
    raw_bam_targets,
    sample_ids,
    transcriptome_bam_targets,
    trimmed_fastq_targets,
    umi_dedup_bam_targets,
)


MUTPROF_TARGETS = []

if module_enabled("fastq_qc"):
    MUTPROF_TARGETS.extend(input_fastqc_targets(config, SAMPLES, RESULTS_DIR))

if module_enabled("trimming"):
    MUTPROF_TARGETS.extend(trimmed_fastq_targets(SAMPLES, RESULTS_DIR, config))

if module_enabled("alignment"):
    MUTPROF_TARGETS.extend(raw_bam_targets(SAMPLES, RESULTS_DIR))

if module_enabled("mark_duplicates"):
    MUTPROF_TARGETS.extend(markdup_bam_targets(SAMPLES, RESULTS_DIR))

if module_enabled("bam_filter") and not module_enabled("umi_deduplication"):
    MUTPROF_TARGETS.extend(filtered_bam_targets(SAMPLES, RESULTS_DIR))

if module_enabled("umi_deduplication"):
    MUTPROF_TARGETS.extend(pre_dedup_bam_targets(SAMPLES, RESULTS_DIR))
    MUTPROF_TARGETS.extend(umi_dedup_bam_targets(SAMPLES, RESULTS_DIR))
    MUTPROF_TARGETS.extend(pre_dedup_qc_targets(SAMPLES, RESULTS_DIR))
    MUTPROF_TARGETS.extend(
        f"{RESULTS_DIR}/qc/umi_tools/{sample}{suffix}"
        for sample in sample_ids(SAMPLES)
        for suffix in (
            "_edit_distance.tsv",
            "_per_umi.tsv",
            "_per_umi_per_position.tsv",
        )
    )

if module_enabled("bam_qc"):
    MUTPROF_TARGETS.extend(bam_qc_targets(SAMPLES, RESULTS_DIR))

if gene_count_mode_enabled("standard"):
    MUTPROF_TARGETS.append(f"{RESULTS_DIR}/counts/featurecounts/gene_counts.txt")

if gene_count_mode_enabled("exon_strict"):
    MUTPROF_TARGETS.append(f"{RESULTS_DIR}/counts/featurecounts/exon_strict_counts.txt")

if gene_count_mode_enabled("full_gene"):
    MUTPROF_TARGETS.append(f"{RESULTS_DIR}/counts/featurecounts/full_gene_counts.txt")

if module_settings("alignment").get("transcriptome_bam", False):
    MUTPROF_TARGETS.extend(transcriptome_bam_targets(SAMPLES, RESULTS_DIR))

if salmon_result_enabled("isoform_results"):
    MUTPROF_TARGETS.extend(
        f"{RESULTS_DIR}/counts/salmon/{sample}/quant.sf" for sample in SAMPLE_IDS
    )

if salmon_result_enabled("gene_results"):
    MUTPROF_TARGETS.extend(
        f"{RESULTS_DIR}/counts/salmon/{sample}/quant.genes.sf" for sample in SAMPLE_IDS
    )

if rsem_result_enabled("isoform_results"):
    MUTPROF_TARGETS.extend(
        f"{RESULTS_DIR}/counts/rsem/{sample}/{sample}.isoforms.results"
        for sample in SAMPLE_IDS
    )

if rsem_result_enabled("gene_results"):
    MUTPROF_TARGETS.extend(
        f"{RESULTS_DIR}/counts/rsem/{sample}/{sample}.genes.results"
        for sample in SAMPLE_IDS
    )

if module_enabled("coverage"):
    MUTPROF_TARGETS.extend(bigwig_targets(SAMPLES, RESULTS_DIR))

if module_enabled("anomaly_calling"):
    MUTPROF_TARGETS.extend(anomaly_targets(RESULTS_DIR))

if module_enabled("mutation_counts"):
    MUTPROF_TARGETS.extend(per_base_targets(SAMPLES, RESULTS_DIR))

if module_enabled("editing_tracks"):
    MUTPROF_TARGETS.extend(editing_track_targets(config, SAMPLES, RESULTS_DIR))

MULTIQC_INPUTS.extend(MUTPROF_TARGETS)

if module_enabled("multiqc"):
    MUTPROF_TARGETS.append(f"{RESULTS_DIR}/multiqc/multiqc_report.html")

ASSAY_TARGETS.extend(MUTPROF_TARGETS)

include: workflow_file("rules/common/rna_quantification.smk")
include: workflow_file("rules/common/mutational_profiling.smk")
