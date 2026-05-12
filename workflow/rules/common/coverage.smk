from params import render_tool_params, resource_value, tool_extra
from tracks import coverage_strand_arg


rule bam_coverage:
    input:
        bam=lambda wildcards: final_bam_path(wildcards.sample),
        bai=lambda wildcards: final_bai_path(wildcards.sample)
    output:
        bw=f"{RESULTS_DIR}/tracks/bigwig/{{sample}}.bw"
    log:
        f"{RESULTS_DIR}/logs/deeptools/bamCoverage.{{sample}}.log"
    threads:
        int(resource_value(config, "deeptools", "threads", 4))
    resources:
        mem_mb=int(resource_value(config, "deeptools", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "deeptools", "runtime_min", 120))
    conda:
        str(WORKFLOW_DIR / "envs" / "deeptools.yaml")
    params:
        rendered=lambda wildcards: render_tool_params(
            config,
            "deeptools",
            overrides={
                "effective_genome_size": config.get("genome", {}).get("effective_genome_size")
            },
        ),
        strand=lambda wildcards: coverage_strand_arg(config, SAMPLES, wildcards.sample),
        extra=lambda wildcards: tool_extra(config, "deeptools")
    shell:
        """
        mkdir -p $(dirname {output.bw}) $(dirname {log})
        bamCoverage -b {input.bam} -o {output.bw} -p {threads} \
            {params.rendered} {params.strand} {params.extra} > {log} 2>&1
        """
