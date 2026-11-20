from clip import crosslink_offset, representative_count_samtools_args, track_strand_args
from params import module_tool_extra, resource_value


rule crosslink_sites:
    input:
        bam=lambda wildcards: final_bam_path(wildcards.sample),
        bai=lambda wildcards: final_bai_path(wildcards.sample)
    output:
        bw=f"{RESULTS_DIR}/tracks/crosslink_sites/{{sample}}.{{strand}}.bw"
    wildcard_constraints:
        strand="plus|minus"
    log:
        f"{RESULTS_DIR}/logs/deeptools/crosslink_sites.{{sample}}.{{strand}}.log"
    threads:
        int(resource_value(config, "deeptools", "threads", 4))
    resources:
        mem_mb=int(resource_value(config, "deeptools", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "deeptools", "runtime_min", 120))
    conda:
        workflow_file("envs/deeptools.yaml")
    params:
        offset=lambda wildcards: crosslink_offset(config),
        strand=lambda wildcards: track_strand_args(
            config, SAMPLES, wildcards.sample, wildcards.strand, "crosslink_sites"
        ),
        min_mapq=lambda wildcards: int(
            module_settings("crosslink_sites").get("min_mapping_quality", 10)
        ),
        extra=lambda wildcards: module_tool_extra(config, "crosslink_sites", "deeptools")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bw}) $(dirname {log})
        bamCoverage -b {input.bam} -o {output.bw} -p {threads} \
            --binSize 1 --normalizeUsing None --Offset {params.offset} \
            --minMappingQuality {params.min_mapq} {params.strand} {params.extra} > {log} 2>&1
        """


rule clip_enrichment:
    input:
        ip_bam=lambda wildcards: final_bam_path(wildcards.sample),
        ip_bai=lambda wildcards: final_bai_path(wildcards.sample),
        input_bam=lambda wildcards: final_bam_path(wildcards.control),
        input_bai=lambda wildcards: final_bai_path(wildcards.control)
    output:
        bw=f"{RESULTS_DIR}/tracks/enrichment/{{sample}}.vs.{{control}}.{{strand}}.bw"
    wildcard_constraints:
        strand="plus|minus"
    log:
        f"{RESULTS_DIR}/logs/deeptools/enrichment.{{sample}}.vs.{{control}}.{{strand}}.log"
    threads:
        int(resource_value(config, "deeptools", "threads", 4))
    resources:
        mem_mb=int(resource_value(config, "deeptools", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "deeptools", "runtime_min", 120))
    conda:
        workflow_file("envs/deeptools.yaml")
    params:
        bin_size=lambda wildcards: int(module_settings("enrichment").get("bin_size", 50)),
        pseudocount=lambda wildcards: float(
            module_settings("enrichment").get("pseudocount", 1)
        ),
        ip_count_args=lambda wildcards: representative_count_samtools_args(
            config, SAMPLES, wildcards.sample, "enrichment"
        ),
        in_count_args=lambda wildcards: representative_count_samtools_args(
            config, SAMPLES, wildcards.control, "enrichment"
        ),
        strand=lambda wildcards: track_strand_args(
            config, SAMPLES, wildcards.sample, wildcards.strand, "enrichment"
        ),
        min_mapq=lambda wildcards: int(
            module_settings("enrichment").get("min_mapping_quality", 10)
        ),
        extra=lambda wildcards: module_tool_extra(config, "enrichment", "deeptools")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bw:q}) $(dirname {log:q})
        ip_reads=$(samtools view -@ {threads} -c {params.ip_count_args} {input.ip_bam:q})
        in_reads=$(samtools view -@ {threads} -c {params.in_count_args} {input.input_bam:q})
        scale_factors=$(python -c 'import sys; n1, n2 = int(sys.argv[1]), int(sys.argv[2]); m = min(n1, n2); print(f"{{m / n1:.12g}}:{{m / n2:.12g}}" if m > 0 else "1:1")' "$ip_reads" "$in_reads")
        bamCompare -b1 {input.ip_bam:q} -b2 {input.input_bam:q} -o {output.bw:q} \
            --numberOfProcessors {threads} --operation log2 \
            --scaleFactorsMethod None --scaleFactors "$scale_factors" \
            --binSize {params.bin_size} \
            --pseudocount {params.pseudocount} {params.pseudocount} \
            --skipZeroOverZero --minMappingQuality {params.min_mapq} \
            {params.strand} {params.extra} > {log:q} 2>&1
        """
