from params import render_tool_params, resource_value, tool_extra
from samples import as_csv


rule bowtie2_align:
    input:
        r1=alignment_r1,
        r2=alignment_r2
    output:
        bam=f"{RESULTS_DIR}/bam/raw/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/bam/raw/{{sample}}.bam.bai"
    log:
        f"{RESULTS_DIR}/logs/bowtie2/{{sample}}.log"
    threads:
        int(resource_value(config, "bowtie2", "threads", 8))
    resources:
        mem_mb=int(resource_value(config, "bowtie2", "mem_mb", 8192)),
        runtime_min=int(resource_value(config, "bowtie2", "runtime_min", 240))
    conda:
        str(WORKFLOW_DIR / "envs" / "bowtie2.yaml")
    params:
        index=lambda wildcards: config["genome"]["bowtie2_index"],
        r1_csv=lambda wildcards, input: as_csv(input.r1),
        r2_csv=lambda wildcards, input: as_csv(input.r2),
        bowtie2_args=lambda wildcards: render_tool_params(config, "bowtie2"),
        sort_args=lambda wildcards: render_tool_params(config, "samtools", section="sort"),
        extra=lambda wildcards: tool_extra(config, "bowtie2")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        bowtie2 -x {params.index} -1 {params.r1_csv} -2 {params.r2_csv} \
            -p {threads} {params.bowtie2_args} {params.extra} 2> {log} \
            | samtools sort -@ {threads} {params.sort_args} -o {output.bam} - 2>> {log}
        samtools index {output.bam} {output.bai} 2>> {log}
        """
