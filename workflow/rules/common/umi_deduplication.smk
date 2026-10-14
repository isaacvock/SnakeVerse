from clip import dedup_args, qname_umi_regex
from params import module_tool_extra, resource_value
from samples import processed_sample_layout


rule tag_qname_umi:
    input:
        bam=lambda wildcards: filtered_bam_path(wildcards.sample)
            if module_enabled("bam_filter")
            else pre_filter_bam_path(wildcards.sample)
    output:
        bam=f"{RESULTS_DIR}/bam/pre_dedup/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/bam/pre_dedup/{{sample}}.bam.bai"
    log:
        f"{RESULTS_DIR}/logs/umi_tools/tag.{{sample}}.log"
    threads: 1
    resources:
        mem_mb=int(resource_value(config, "umi_tools_tag", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "umi_tools_tag", "runtime_min", 120))
    conda:
        workflow_file("envs/umi_tools.yaml")
    params:
        regex=lambda wildcards: qname_umi_regex(config),
        script=lambda wildcards: workflow_file("scripts/qname_umi_to_rx.py")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        python {params.script} \
            --input {input.bam} --output {output.bam} --regex {params.regex:q} > {log} 2>&1
        samtools index {output.bam} {output.bai} >> {log} 2>&1
        """


rule umi_tools_dedup:
    input:
        bam=f"{RESULTS_DIR}/bam/pre_dedup/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/bam/pre_dedup/{{sample}}.bam.bai"
    output:
        bam=f"{RESULTS_DIR}/bam/umi_dedup/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/bam/umi_dedup/{{sample}}.bam.bai",
        edit_distance=f"{RESULTS_DIR}/qc/umi_tools/{{sample}}_edit_distance.tsv",
        per_umi=f"{RESULTS_DIR}/qc/umi_tools/{{sample}}_per_umi.tsv",
        per_position=f"{RESULTS_DIR}/qc/umi_tools/{{sample}}_per_umi_per_position.tsv"
    log:
        f"{RESULTS_DIR}/logs/umi_tools/dedup.{{sample}}.log"
    threads: 1
    resources:
        mem_mb=int(resource_value(config, "umi_tools_dedup", "mem_mb", 8192)),
        runtime_min=int(resource_value(config, "umi_tools_dedup", "runtime_min", 240))
    conda:
        workflow_file("envs/umi_tools.yaml")
    params:
        rendered=lambda wildcards: dedup_args(
            config, processed_sample_layout(config, SAMPLES, wildcards.sample) == "paired"
        ),
        stats=lambda wildcards: f"{RESULTS_DIR}/qc/umi_tools/{wildcards.sample}",
        extra=lambda wildcards: module_tool_extra(config, "umi_deduplication", "umi_tools")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {output.edit_distance}) $(dirname {log})
        umi_tools dedup --stdin={input.bam} --stdout={output.bam} \
            {params.rendered} {params.extra} --output-stats={params.stats} --log={log}
        samtools index {output.bam} {output.bai} >> {log} 2>&1
        """


rule samtools_prededup_qc:
    input:
        bam=f"{RESULTS_DIR}/bam/pre_dedup/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/bam/pre_dedup/{{sample}}.bam.bai"
    output:
        flagstat=f"{RESULTS_DIR}/qc/bam/pre_dedup/{{sample}}.flagstat.txt",
        idxstats=f"{RESULTS_DIR}/qc/bam/pre_dedup/{{sample}}.idxstats.txt"
    wildcard_constraints:
        sample="[^/]+"
    log:
        f"{RESULTS_DIR}/logs/samtools/qc.pre_dedup.{{sample}}.log"
    threads:
        int(resource_value(config, "samtools_qc", "threads", 2))
    resources:
        mem_mb=int(resource_value(config, "samtools_qc", "mem_mb", 2048)),
        runtime_min=int(resource_value(config, "samtools_qc", "runtime_min", 60))
    conda:
        workflow_file("envs/samtools.yaml")
    shell:
        """
        mkdir -p $(dirname {output.flagstat}) $(dirname {log})
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat} 2> {log}
        samtools idxstats {input.bam} > {output.idxstats} 2>> {log}
        """
