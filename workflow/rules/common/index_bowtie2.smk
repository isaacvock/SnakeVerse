from params import render_tool_params, resource_value, tool_extra
from refs import generated_index_prefix, genome_fasta


rule bowtie2_index:
    input:
        fasta=lambda wildcards: genome_fasta(config)
    output:
        touch(f"{RESULTS_DIR}/reference/bowtie2/{GENOME_SLUG}/.snakeverse_bowtie2_index.done")
    log:
        f"{RESULTS_DIR}/logs/bowtie2/build_index.log"
    threads:
        int(resource_value(config, "bowtie2_index", "threads", 2))
    resources:
        mem_mb=int(resource_value(config, "bowtie2_index", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "bowtie2_index", "runtime_min", 120))
    conda:
        str(WORKFLOW_DIR / "envs" / "bowtie2.yaml")
    params:
        prefix=lambda wildcards: generated_index_prefix(config, RESULTS_DIR, "bowtie2"),
        rendered=lambda wildcards: render_tool_params(config, "bowtie2", section="index"),
        extra=lambda wildcards: tool_extra(config, "bowtie2")
    shell:
        """
        mkdir -p $(dirname {params.prefix}) $(dirname {log})
        bowtie2-build --threads {threads} {params.rendered} {params.extra} \
            {input.fasta} {params.prefix} > {log} 2>&1
        """

