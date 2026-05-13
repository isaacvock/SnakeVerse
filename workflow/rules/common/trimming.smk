from params import render_tool_params, resource_value, tool_extra
from samples import fastq_for_read, optional_fastq_for_read, unit_layout


if TRIMMING_TOOL == "fastp":
    rule fastp:
        input:
            r1=lambda wildcards: fastq_for_read(SAMPLES, wildcards.unit, "R1", RESULTS_DIR),
            r2=lambda wildcards: optional_fastq_for_read(SAMPLES, wildcards.unit, "R2", RESULTS_DIR)
        output:
            r1=f"{RESULTS_DIR}/fastq/trimmed/{{unit}}_R1.fastq.gz",
            r2=f"{RESULTS_DIR}/fastq/trimmed/{{unit}}_R2.fastq.gz",
            html=f"{RESULTS_DIR}/qc/fastp/{{unit}}.html",
            json=f"{RESULTS_DIR}/qc/fastp/{{unit}}.json"
        log:
            f"{RESULTS_DIR}/logs/fastp/{{unit}}.log"
        threads:
            int(resource_value(config, "fastp", "threads", 4))
        resources:
            mem_mb=int(resource_value(config, "fastp", "mem_mb", 4096)),
            runtime_min=int(resource_value(config, "fastp", "runtime_min", 120))
        conda:
            str(WORKFLOW_DIR / "envs" / "fastp.yaml")
        params:
            layout=lambda wildcards: unit_layout(SAMPLES, wildcards.unit),
            rendered=lambda wildcards: render_tool_params(config, "fastp"),
            extra=lambda wildcards: tool_extra(config, "fastp")
        shell:
            """
            mkdir -p $(dirname {output.r1}) $(dirname {output.html}) $(dirname {log})
            if [ "{params.layout}" = "paired" ]; then
                fastp {params.rendered} {params.extra} --thread {threads} \
                    --in1 {input.r1} --in2 {input.r2} \
                    --out1 {output.r1} --out2 {output.r2} \
                    --html {output.html} --json {output.json} > {log} 2>&1
            else
                fastp {params.rendered} {params.extra} --thread {threads} \
                    --in1 {input.r1} --out1 {output.r1} \
                    --html {output.html} --json {output.json} > {log} 2>&1
                : > {output.r2}
            fi
            """
elif TRIMMING_TOOL == "cutadapt":
    rule cutadapt:
        input:
            r1=lambda wildcards: fastq_for_read(SAMPLES, wildcards.unit, "R1", RESULTS_DIR),
            r2=lambda wildcards: optional_fastq_for_read(SAMPLES, wildcards.unit, "R2", RESULTS_DIR)
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
else:
    raise ValueError(f"Unsupported trimming.tool: {TRIMMING_TOOL}")
