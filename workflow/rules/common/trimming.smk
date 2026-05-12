from params import render_tool_params, resource_value, tool_extra
from samples import fastq_for_read, optional_fastq_for_read, unit_layout


rule cutadapt:
    input:
        r1=lambda wildcards: fastq_for_read(SAMPLES, wildcards.unit, "R1"),
        r2=lambda wildcards: optional_fastq_for_read(SAMPLES, wildcards.unit, "R2")
    output:
        r1=f"{RESULTS_DIR}/fastq/trimmed/{{unit}}_R1.fastq.gz",
        r2=f"{RESULTS_DIR}/fastq/trimmed/{{unit}}_R2.fastq.gz"
    log:
        f"{RESULTS_DIR}/logs/cutadapt/{{unit}}.log"
    threads:
        int(resource_value(config, "cutadapt", "threads", 4))
    resources:
        mem_mb=int(resource_value(config, "cutadapt", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "cutadapt", "runtime_min", 120))
    conda:
        str(WORKFLOW_DIR / "envs" / "cutadapt.yaml")
    params:
        layout=lambda wildcards: unit_layout(SAMPLES, wildcards.unit),
        rendered=lambda wildcards: render_tool_params(config, "cutadapt"),
        extra=lambda wildcards: tool_extra(config, "cutadapt")
    shell:
        """
        mkdir -p $(dirname {output.r1}) $(dirname {log})
        if [ "{params.layout}" = "paired" ]; then
            cutadapt {params.rendered} {params.extra} --cores {threads} \
                -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log} 2>&1
        else
            cutadapt {params.rendered} {params.extra} --cores {threads} \
                -o {output.r1} {input.r1} > {log} 2>&1
            : > {output.r2}
        fi
        """
