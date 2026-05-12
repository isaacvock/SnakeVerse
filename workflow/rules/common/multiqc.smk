from params import render_tool_params, resource_value, tool_extra


rule multiqc:
    input:
        lambda wildcards: MULTIQC_INPUTS
    output:
        html=f"{RESULTS_DIR}/multiqc/multiqc_report.html"
    log:
        f"{RESULTS_DIR}/logs/multiqc.log"
    threads:
        int(resource_value(config, "multiqc", "threads", 1))
    resources:
        mem_mb=int(resource_value(config, "multiqc", "mem_mb", 2048)),
        runtime_min=int(resource_value(config, "multiqc", "runtime_min", 60))
    conda:
        str(WORKFLOW_DIR / "envs" / "multiqc.yaml")
    params:
        search_dir=lambda wildcards, output: Path(output.html).parents[1].as_posix(),
        outdir=lambda wildcards, output: Path(output.html).parent.as_posix(),
        filename=lambda wildcards, output: Path(output.html).name,
        rendered=lambda wildcards: render_tool_params(config, "multiqc"),
        extra=lambda wildcards: tool_extra(config, "multiqc")
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        multiqc {params.search_dir} --outdir {params.outdir} --filename {params.filename} \
            {params.rendered} {params.extra} > {log} 2>&1
        """
