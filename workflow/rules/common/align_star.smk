from params import render_tool_params, resource_value, tool_extra
from samples import as_csv


rule star_align:
    input:
        r1=alignment_r1,
        r2=alignment_r2
    output:
        bam=f"{RESULTS_DIR}/bam/raw/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/bam/raw/{{sample}}.bam.bai"
    log:
        f"{RESULTS_DIR}/logs/star/{{sample}}.log"
    threads:
        int(resource_value(config, "star", "threads", 12))
    resources:
        mem_mb=int(resource_value(config, "star", "mem_mb", 32000)),
        runtime_min=int(resource_value(config, "star", "runtime_min", 360))
    conda:
        str(WORKFLOW_DIR / "envs" / "star.yaml")
    params:
        genome_dir=lambda wildcards: config["genome"]["star_index"],
        r1_csv=lambda wildcards, input: as_csv(input.r1),
        r2_csv=lambda wildcards, input: as_csv(input.r2),
        prefix=lambda wildcards: f"{RESULTS_DIR}/star/{wildcards.sample}/",
        star_args=lambda wildcards: render_tool_params(config, "star"),
        extra=lambda wildcards: tool_extra(config, "star")
    shell:
        """
        set -euo pipefail
        mkdir -p {params.prefix} $(dirname {output.bam}) $(dirname {log})
        STAR --runThreadN {threads} \
            --genomeDir {params.genome_dir} \
            --readFilesIn {params.r1_csv} {params.r2_csv} \
            --outFileNamePrefix {params.prefix} \
            {params.star_args} {params.extra} > {log} 2>&1
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
        samtools index {output.bam} {output.bai} 2>> {log}
        """
