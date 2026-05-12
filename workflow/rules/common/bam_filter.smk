from params import render_tool_params, resource_value, tool_extra


rule samtools_filter:
    input:
        bam=lambda wildcards: raw_bam_path(wildcards.sample),
        bai=lambda wildcards: raw_bai_path(wildcards.sample)
    output:
        bam=f"{RESULTS_DIR}/bam/filtered/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/bam/filtered/{{sample}}.bam.bai"
    log:
        f"{RESULTS_DIR}/logs/samtools/filter.{{sample}}.log"
    threads:
        int(resource_value(config, "samtools_filter", "threads", 4))
    resources:
        mem_mb=int(resource_value(config, "samtools_filter", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "samtools_filter", "runtime_min", 120))
    conda:
        str(WORKFLOW_DIR / "envs" / "samtools.yaml")
    params:
        view_args=lambda wildcards: render_tool_params(config, "samtools", section="filter"),
        extra=lambda wildcards: tool_extra(config, "samtools")
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        samtools view -@ {threads} -b {params.view_args} {params.extra} \
            -o {output.bam} {input.bam} > {log} 2>&1
        samtools index {output.bam} {output.bai} >> {log} 2>&1
        """
