from params import resource_value
from samples import sample_layout


rule samtools_mark_duplicates:
    input:
        bam=lambda wildcards: raw_bam_path(wildcards.sample),
        bai=lambda wildcards: raw_bai_path(wildcards.sample)
    output:
        bam=f"{RESULTS_DIR}/bam/markdup/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/bam/markdup/{{sample}}.bam.bai"
    log:
        f"{RESULTS_DIR}/logs/samtools/markdup.{{sample}}.log"
    threads:
        int(resource_value(config, "samtools_markdup", "threads", 4))
    resources:
        mem_mb=int(resource_value(config, "samtools_markdup", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "samtools_markdup", "runtime_min", 120))
    conda:
        str(WORKFLOW_DIR / "envs" / "samtools.yaml")
    params:
        layout=lambda wildcards: sample_layout(SAMPLES, wildcards.sample)
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        tmpdir=$(mktemp -d)
        if [ "{params.layout}" = "paired" ]; then
            samtools sort -@ {threads} -n -o "$tmpdir/{wildcards.sample}.name.bam" {input.bam} > {log} 2>&1
            samtools fixmate -@ {threads} -m "$tmpdir/{wildcards.sample}.name.bam" "$tmpdir/{wildcards.sample}.fixmate.bam" >> {log} 2>&1
            samtools sort -@ {threads} -o "$tmpdir/{wildcards.sample}.coord.bam" "$tmpdir/{wildcards.sample}.fixmate.bam" >> {log} 2>&1
            samtools markdup -@ {threads} "$tmpdir/{wildcards.sample}.coord.bam" {output.bam} >> {log} 2>&1
        else
            samtools markdup -@ {threads} {input.bam} {output.bam} > {log} 2>&1
        fi
        samtools index {output.bam} {output.bai} >> {log} 2>&1
        rm -rf "$tmpdir"
        """
