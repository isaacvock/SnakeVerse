from params import render_tool_params, resource_value, tool_extra
from refs import aligner_index_inputs, aligner_index_prefix
from samples import bowtie2_reads_arg


rule bowtie2_align:
    input:
        index=lambda wildcards: aligner_index_inputs(config, RESULTS_DIR, "bowtie2"),
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
        index=lambda wildcards: aligner_index_prefix(config, RESULTS_DIR, "bowtie2"),
        reads=lambda wildcards, input: bowtie2_reads_arg(input.r1, input.r2),
        bowtie2_args=lambda wildcards: render_tool_params(config, "bowtie2", section="align"),
        sort_args=lambda wildcards: render_tool_params(config, "samtools", section="sort"),
        extra=lambda wildcards: tool_extra(config, "bowtie2")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        bowtie2 -x {params.index} {params.reads} \
            -p {threads} {params.bowtie2_args} {params.extra} 2> {log} \
            | samtools sort -@ {threads} {params.sort_args} -o {output.bam} - 2>> {log}
        samtools index {output.bam} {output.bai} 2>> {log}
        """
