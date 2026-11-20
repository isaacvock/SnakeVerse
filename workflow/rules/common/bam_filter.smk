from params import render_tool_params, resource_value, tool_extra
from samples import processed_sample_layout


rule samtools_filter:
    input:
        bam=lambda wildcards: pre_filter_bam_path(wildcards.sample),
        bai=lambda wildcards: pre_filter_bai_path(wildcards.sample)
    output:
        bam=temp(f"{RESULTS_DIR}/bam/filtered/{{sample}}.bam")
            if module_enabled("umi_deduplication")
            else f"{RESULTS_DIR}/bam/filtered/{{sample}}.bam",
        bai=temp(f"{RESULTS_DIR}/bam/filtered/{{sample}}.bam.bai")
            if module_enabled("umi_deduplication")
            else f"{RESULTS_DIR}/bam/filtered/{{sample}}.bam.bai"
    log:
        f"{RESULTS_DIR}/logs/samtools/filter.{{sample}}.log"
    threads:
        int(resource_value(config, "samtools_filter", "threads", 4))
    resources:
        mem_mb=int(resource_value(config, "samtools_filter", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "samtools_filter", "runtime_min", 120))
    conda:
        workflow_file("envs/samtools.yaml")
    params:
        view_args=lambda wildcards: render_tool_params(
            config,
            "samtools",
            section="filter",
            overrides={
                "min_mapq": config.get("modules", {}).get("bam_filter", {}).get("min_mapq"),
                "required_flags": config.get("modules", {}).get("bam_filter", {}).get("required_flags"),
                "excluded_flags": config.get("modules", {}).get("bam_filter", {}).get("excluded_flags"),
                "keep_duplicates": config.get("modules", {}).get("bam_filter", {}).get("keep_duplicates"),
            },
        ),
        excluded_contigs=lambda wildcards: ",".join(config.get("modules", {}).get("bam_filter", {}).get("exclude_contigs", []) or []),
        blacklist=lambda wildcards: config.get("reference", {}).get("blacklist", "") or "",
        paired=lambda wildcards: processed_sample_layout(config, SAMPLES, wildcards.sample) == "paired",
        extra=lambda wildcards: tool_extra(config, "samtools", section="filter")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bam:q}) $(dirname {log:q})
        tmp_bam=$(mktemp --tmpdir=$(dirname {output.bam:q}) {wildcards.sample}.filter.XXXXXX.bam)
        bad_qnames="$tmp_bam.bad_qnames"
        trap 'rm -f "$tmp_bam" "$bad_qnames"' EXIT
        samtools view -@ {threads} -h {params.view_args} {params.extra} {input.bam:q} \
            | awk -v exclude="{params.excluded_contigs}" 'BEGIN {{
                n = split(exclude, contigs, ",");
                for (i = 1; i <= n; i++) if (contigs[i] != "") skip[contigs[i]] = 1;
            }}
            /^@/ {{ print; next }}
            !($3 in skip) {{ print }}' \
            | samtools view -@ {threads} -b -o "$tmp_bam" - > {log:q} 2>&1
        if [ -n "{params.blacklist}" ]; then
            if [ "{params.paired}" = "True" ]; then
                bedtools intersect -u -abam "$tmp_bam" -b {params.blacklist:q} 2>> {log:q} \
                    | samtools view - \
                    | cut -f1 \
                    | LC_ALL=C sort -u > "$bad_qnames"
                if [ -s "$bad_qnames" ]; then
                    samtools view -@ {threads} -b -N "^$bad_qnames" -o {output.bam:q} "$tmp_bam" >> {log:q} 2>&1
                else
                    mv "$tmp_bam" {output.bam:q}
                fi
            else
                bedtools intersect -v -abam "$tmp_bam" -b {params.blacklist:q} > {output.bam:q} 2>> {log:q}
            fi
        else
            mv "$tmp_bam" {output.bam:q}
        fi
        samtools index {output.bam:q} {output.bai:q} >> {log:q} 2>&1
        """
