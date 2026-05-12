GENERIC_TARGETS = []

if output_enabled("fastq_qc"):
    GENERIC_TARGETS.extend(fastqc_targets(SAMPLES, RESULTS_DIR))

if output_enabled("trimmed_fastq") and step_enabled("trimming", False):
    GENERIC_TARGETS.extend(trimmed_fastq_targets(SAMPLES, RESULTS_DIR))

if output_enabled("raw_bam"):
    GENERIC_TARGETS.extend(raw_bam_targets(SAMPLES, RESULTS_DIR))

if output_enabled("filtered_bam") and step_enabled("bam_filter", True):
    GENERIC_TARGETS.extend(filtered_bam_targets(SAMPLES, RESULTS_DIR))

GENERIC_TARGETS.extend(bam_qc_targets(SAMPLES, RESULTS_DIR))

if output_enabled("bigwig") and step_enabled("coverage", False):
    GENERIC_TARGETS.extend(bigwig_targets(SAMPLES, RESULTS_DIR))

MULTIQC_INPUTS.extend(GENERIC_TARGETS)

if output_enabled("multiqc"):
    GENERIC_TARGETS.append(f"{RESULTS_DIR}/multiqc/multiqc_report.html")

ASSAY_TARGETS.extend(GENERIC_TARGETS)

