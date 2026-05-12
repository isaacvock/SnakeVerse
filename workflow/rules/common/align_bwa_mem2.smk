from params import render_tool_params, resource_value, tool_extra
from refs import aligner_index_inputs, aligner_index_prefix
from samples import bwa_mem2_reads_arg


rule bwa_mem2_align:
    input:
        index=lambda wildcards: aligner_index_inputs(config, RESULTS_DIR, "bwa_mem2"),
        r1=alignment_r1,
        r2=alignment_r2
    output:
        bam=f"{RESULTS_DIR}/bam/raw/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/bam/raw/{{sample}}.bam.bai"
    log:
        f"{RESULTS_DIR}/logs/bwa_mem2/{{sample}}.log"
    threads:
        int(resource_value(config, "bwa_mem2", "threads", 8))
    resources:
        mem_mb=int(resource_value(config, "bwa_mem2", "mem_mb", 8192)),
        runtime_min=int(resource_value(config, "bwa_mem2", "runtime_min", 240))
    conda:
        str(WORKFLOW_DIR / "envs" / "bwa_mem2.yaml")
    params:
        index=lambda wildcards: aligner_index_prefix(config, RESULTS_DIR, "bwa_mem2"),
        reads=lambda wildcards, input: bwa_mem2_reads_arg(input.r1, input.r2),
        bwa_args=lambda wildcards: render_tool_params(config, "bwa_mem2", section="align"),
        sort_args=lambda wildcards: render_tool_params(config, "samtools", section="sort"),
        extra=lambda wildcards: tool_extra(config, "bwa_mem2")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        bwa-mem2 mem -t {threads} {params.bwa_args} {params.extra} \
            {params.index} {params.reads} 2> {log} \
            | samtools sort -@ {threads} {params.sort_args} -o {output.bam} - 2>> {log}
        samtools index {output.bam} {output.bai} 2>> {log}
        """

