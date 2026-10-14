from clip import crosslink_offset, track_strand_args
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
        strand=lambda wildcards: track_strand_args(
            config, SAMPLES, wildcards.sample, wildcards.strand, "enrichment"
        ),
        min_mapq=lambda wildcards: int(
            module_settings("enrichment").get("min_mapping_quality", 10)
        ),
        extra=lambda wildcards: module_tool_extra(config, "enrichment", "deeptools")
    shell:
        """
        mkdir -p $(dirname {output.bw}) $(dirname {log})
        bamCompare -b1 {input.ip_bam} -b2 {input.input_bam} -o {output.bw} \
            --numberOfProcessors {threads} --operation log2 \
            --scaleFactorsMethod readCount --binSize {params.bin_size} \
            --pseudocount {params.pseudocount} {params.pseudocount} \
            --skipZeroOverZero --minMappingQuality {params.min_mapq} \
            {params.strand} {params.extra} > {log} 2>&1
        """
