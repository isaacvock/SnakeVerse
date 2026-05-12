from params import render_tool_params, resource_value, tool_extra
from refs import aligner_index_inputs, aligner_index_prefix
from samples import star_reads_arg


if output_enabled("transcriptome_bam"):
    rule star_align:
        input:
            index=lambda wildcards: aligner_index_inputs(config, RESULTS_DIR, "star"),
            r1=alignment_r1,
            r2=alignment_r2
        output:
            bam=f"{RESULTS_DIR}/bam/raw/{{sample}}.bam",
            bai=f"{RESULTS_DIR}/bam/raw/{{sample}}.bam.bai",
            transcriptome=f"{RESULTS_DIR}/bam/transcriptome/{{sample}}.bam"
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
            genome_dir=lambda wildcards: aligner_index_prefix(config, RESULTS_DIR, "star"),
            reads=lambda wildcards, input: star_reads_arg(input.r1, input.r2),
            prefix=lambda wildcards: f"{RESULTS_DIR}/star/{wildcards.sample}/",
            star_args=lambda wildcards: render_tool_params(config, "star", section="align"),
            extra=lambda wildcards: tool_extra(config, "star")
        shell:
            """
            set -euo pipefail
            mkdir -p {params.prefix} $(dirname {output.bam}) $(dirname {output.transcriptome}) $(dirname {log})
            STAR --runThreadN {threads} \
                --genomeDir {params.genome_dir} \
                --readFilesIn {params.reads} \
                --outFileNamePrefix {params.prefix} \
                {params.star_args} {params.extra} > {log} 2>&1
            mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
            mv {params.prefix}Aligned.toTranscriptome.out.bam {output.transcriptome}
            samtools index {output.bam} {output.bai} 2>> {log}
            """
else:
    rule star_align:
        input:
            index=lambda wildcards: aligner_index_inputs(config, RESULTS_DIR, "star"),
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
            genome_dir=lambda wildcards: aligner_index_prefix(config, RESULTS_DIR, "star"),
            reads=lambda wildcards, input: star_reads_arg(input.r1, input.r2),
            prefix=lambda wildcards: f"{RESULTS_DIR}/star/{wildcards.sample}/",
            star_args=lambda wildcards: render_tool_params(config, "star", section="align"),
            extra=lambda wildcards: tool_extra(config, "star")
        shell:
            """
            set -euo pipefail
            mkdir -p {params.prefix} $(dirname {output.bam}) $(dirname {log})
            STAR --runThreadN {threads} \
                --genomeDir {params.genome_dir} \
                --readFilesIn {params.reads} \
                --outFileNamePrefix {params.prefix} \
                {params.star_args} {params.extra} > {log} 2>&1
            mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
            samtools index {output.bam} {output.bai} 2>> {log}
            """
