from params import render_tool_params, resource_value, tool_extra
from refs import generated_index_prefix, genome_fasta


rule bwa_mem2_index:
    input:
        fasta=lambda wildcards: genome_fasta(config)
    output:
        touch(f"{RESULTS_DIR}/reference/bwa_mem2/{GENOME_SLUG}/.snakeverse_bwa_mem2_index.done")
    log:
        f"{RESULTS_DIR}/logs/bwa_mem2/build_index.log"
    threads:
        int(resource_value(config, "bwa_mem2_index", "threads", 2))
    resources:
        mem_mb=int(resource_value(config, "bwa_mem2_index", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "bwa_mem2_index", "runtime_min", 120))
    conda:
        str(WORKFLOW_DIR / "envs" / "bwa_mem2.yaml")
    params:
        prefix=lambda wildcards: generated_index_prefix(config, RESULTS_DIR, "bwa_mem2"),
        rendered=lambda wildcards: render_tool_params(config, "bwa_mem2", section="index"),
        extra=lambda wildcards: tool_extra(config, "bwa_mem2")
    shell:
        """
        mkdir -p $(dirname {params.prefix}) $(dirname {log})
        bwa-mem2 index {params.rendered} {params.extra} -p {params.prefix} \
            {input.fasta} > {log} 2>&1
        """

