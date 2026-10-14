from samples import (
    bam_qc_targets,
    clip_coverage_targets,
    crosslink_site_targets,
    enrichment_targets,
    input_fastqc_targets,
    pool_ids,
    pre_dedup_bam_targets,
    pre_dedup_qc_targets,
    raw_bam_targets,
    sample_ids,
    trimmed_fastq_targets,
    umi_dedup_bam_targets,
)


ECLIP_TARGETS = []

if module_enabled("fastq_qc"):
    ECLIP_TARGETS.extend(input_fastqc_targets(config, SAMPLES, RESULTS_DIR))

if demultiplexing_active(config):
    ECLIP_TARGETS.extend(
        f"{RESULTS_DIR}/qc/cutadapt/demultiplex.{pool}.json"
        for pool in pool_ids(SAMPLES, barcoded_only=True)
    )

if module_enabled("trimming"):
    ECLIP_TARGETS.extend(trimmed_fastq_targets(SAMPLES, RESULTS_DIR, config))

if module_enabled("alignment"):
    ECLIP_TARGETS.extend(raw_bam_targets(SAMPLES, RESULTS_DIR))

if module_enabled("umi_deduplication"):
    ECLIP_TARGETS.extend(pre_dedup_bam_targets(SAMPLES, RESULTS_DIR))
    ECLIP_TARGETS.extend(umi_dedup_bam_targets(SAMPLES, RESULTS_DIR))
    ECLIP_TARGETS.extend(pre_dedup_qc_targets(SAMPLES, RESULTS_DIR))
    ECLIP_TARGETS.extend(
        f"{RESULTS_DIR}/qc/umi_tools/{sample}{suffix}"
        for sample in sample_ids(SAMPLES)
        for suffix in (
            "_edit_distance.tsv",
            "_per_umi.tsv",
            "_per_umi_per_position.tsv",
        )
    )

if module_enabled("bam_qc"):
    ECLIP_TARGETS.extend(bam_qc_targets(SAMPLES, RESULTS_DIR))

if module_enabled("coverage"):
    ECLIP_TARGETS.extend(clip_coverage_targets(SAMPLES, RESULTS_DIR))

if module_enabled("crosslink_sites"):
    ECLIP_TARGETS.extend(crosslink_site_targets(SAMPLES, RESULTS_DIR))

if module_enabled("enrichment"):
    ECLIP_TARGETS.extend(enrichment_targets(SAMPLES, RESULTS_DIR))

MULTIQC_INPUTS.extend(ECLIP_TARGETS)

if module_enabled("multiqc"):
    ECLIP_TARGETS.append(f"{RESULTS_DIR}/multiqc/multiqc_report.html")

ASSAY_TARGETS.extend(ECLIP_TARGETS)
