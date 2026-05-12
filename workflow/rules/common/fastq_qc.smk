from params import render_tool_params, resource_value, tool_extra
from samples import fastq_for_read


rule fastqc:
    input:
        lambda wildcards: fastq_for_read(SAMPLES, wildcards.unit, wildcards.read)
    output:
        directory(f"{RESULTS_DIR}/qc/fastqc/{{unit}}.{{read}}")
    log:
        f"{RESULTS_DIR}/logs/fastqc/{{unit}}.{{read}}.log"
    threads:
        int(resource_value(config, "fastqc", "threads", 1))
    resources:
        mem_mb=int(resource_value(config, "fastqc", "mem_mb", 1024)),
        runtime_min=int(resource_value(config, "fastqc", "runtime_min", 30))
    conda:
        str(WORKFLOW_DIR / "envs" / "fastqc.yaml")
    params:
        rendered=lambda wildcards: render_tool_params(config, "fastqc"),
        extra=lambda wildcards: tool_extra(config, "fastqc")
    shell:
        """
        mkdir -p {output} $(dirname {log})
        fastqc --threads {threads} {params.rendered} {params.extra} --outdir {output} {input} > {log} 2>&1
        """
