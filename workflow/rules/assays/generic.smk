GENERIC_TARGETS = []

if module_enabled("fastq_qc"):
    GENERIC_TARGETS.extend(fastqc_targets(SAMPLES, RESULTS_DIR))

if module_enabled("trimming"):
    GENERIC_TARGETS.extend(trimmed_fastq_targets(SAMPLES, RESULTS_DIR))

if module_enabled("alignment"):
    GENERIC_TARGETS.extend(raw_bam_targets(SAMPLES, RESULTS_DIR))

if module_enabled("mark_duplicates"):
    GENERIC_TARGETS.extend(markdup_bam_targets(SAMPLES, RESULTS_DIR))

if module_enabled("bam_filter"):
    GENERIC_TARGETS.extend(filtered_bam_targets(SAMPLES, RESULTS_DIR))

if module_enabled("bam_qc"):
    GENERIC_TARGETS.extend(bam_qc_targets(SAMPLES, RESULTS_DIR))

if module_enabled("coverage"):
    GENERIC_TARGETS.extend(bigwig_targets(SAMPLES, RESULTS_DIR))

MULTIQC_INPUTS.extend(GENERIC_TARGETS)

if module_enabled("multiqc"):
    GENERIC_TARGETS.append(f"{RESULTS_DIR}/multiqc/multiqc_report.html")

ASSAY_TARGETS.extend(GENERIC_TARGETS)
