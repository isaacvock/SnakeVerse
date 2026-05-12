from params import render_tool_params, resource_value, tool_extra
from refs import generated_index_dir, genome_fasta, star_gtf_arg


rule star_index:
    input:
        fasta=lambda wildcards: genome_fasta(config)
    output:
        touch(f"{RESULTS_DIR}/reference/star/{GENOME_SLUG}/.snakeverse_star_index.done")
    log:
        f"{RESULTS_DIR}/logs/star/build_index.log"
    threads:
        int(resource_value(config, "star_index", "threads", 4))
    resources:
        mem_mb=int(resource_value(config, "star_index", "mem_mb", 16000)),
        runtime_min=int(resource_value(config, "star_index", "runtime_min", 240))
    conda:
        str(WORKFLOW_DIR / "envs" / "star.yaml")
    params:
        genome_dir=lambda wildcards: generated_index_dir(config, RESULTS_DIR, "star"),
        gtf_arg=lambda wildcards: star_gtf_arg(config),
        rendered=lambda wildcards: render_tool_params(config, "star", section="index"),
        extra=lambda wildcards: tool_extra(config, "star")
    shell:
        """
        mkdir -p {params.genome_dir} $(dirname {log})
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {params.genome_dir} \
            --genomeFastaFiles {input.fasta} \
            {params.gtf_arg} {params.rendered} {params.extra} > {log} 2>&1
        """

