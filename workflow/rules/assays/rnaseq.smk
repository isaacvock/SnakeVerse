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

include: workflow_file("rules/common/rna_quantification.smk")
