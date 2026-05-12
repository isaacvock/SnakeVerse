from params import resource_value


rule samtools_bam_qc:
    input:
        bam=lambda wildcards: final_bam_path(wildcards.sample),
        bai=lambda wildcards: final_bai_path(wildcards.sample)
    output:
        flagstat=f"{RESULTS_DIR}/qc/bam/{{sample}}.flagstat.txt",
        idxstats=f"{RESULTS_DIR}/qc/bam/{{sample}}.idxstats.txt"
    log:
        f"{RESULTS_DIR}/logs/samtools/qc.{{sample}}.log"
    threads:
        int(resource_value(config, "samtools_qc", "threads", 2))
    resources:
        mem_mb=int(resource_value(config, "samtools_qc", "mem_mb", 2048)),
        runtime_min=int(resource_value(config, "samtools_qc", "runtime_min", 60))
    conda:
        str(WORKFLOW_DIR / "envs" / "samtools.yaml")
    shell:
        """
        mkdir -p $(dirname {output.flagstat}) $(dirname {log})
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat} 2> {log}
        samtools idxstats {input.bam} > {output.idxstats} 2>> {log}
        """
