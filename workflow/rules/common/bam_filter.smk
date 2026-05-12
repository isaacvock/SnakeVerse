from params import render_tool_params, resource_value, tool_extra


rule samtools_filter:
    input:
        bam=lambda wildcards: pre_filter_bam_path(wildcards.sample),
        bai=lambda wildcards: pre_filter_bai_path(wildcards.sample)
    output:
        bam=f"{RESULTS_DIR}/bam/filtered/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/bam/filtered/{{sample}}.bam.bai"
    log:
        f"{RESULTS_DIR}/logs/samtools/filter.{{sample}}.log"
    threads:
        int(resource_value(config, "samtools_filter", "threads", 4))
    resources:
        mem_mb=int(resource_value(config, "samtools_filter", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "samtools_filter", "runtime_min", 120))
    conda:
        str(WORKFLOW_DIR / "envs" / "samtools.yaml")
    params:
        view_args=lambda wildcards: render_tool_params(
            config,
            "samtools",
            section="filter",
            overrides={
                "min_mapq": config.get("filtering", {}).get("min_mapq"),
                "required_flags": config.get("filtering", {}).get("required_flags"),
                "excluded_flags": config.get("filtering", {}).get("excluded_flags"),
                "keep_duplicates": config.get("filtering", {}).get("keep_duplicates"),
            },
        ),
        excluded_contigs=lambda wildcards: ",".join(config.get("filtering", {}).get("exclude_contigs", []) or []),
        blacklist=lambda wildcards: config.get("genome", {}).get("blacklist", "") or "",
        extra=lambda wildcards: tool_extra(config, "samtools")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        tmp_bam=$(mktemp --tmpdir=$(dirname {output.bam}) {wildcards.sample}.filter.XXXXXX.bam)
        samtools view -@ {threads} -h {params.view_args} {params.extra} {input.bam} \
            | awk -v exclude="{params.excluded_contigs}" 'BEGIN {{
                n = split(exclude, contigs, ",");
                for (i = 1; i <= n; i++) if (contigs[i] != "") skip[contigs[i]] = 1;
            }}
            /^@/ {{ print; next }}
            !($3 in skip) {{ print }}' \
            | samtools view -@ {threads} -b -o "$tmp_bam" - > {log} 2>&1
        if [ -n "{params.blacklist}" ]; then
            bedtools intersect -v -abam "$tmp_bam" -b {params.blacklist} > {output.bam} 2>> {log}
            rm -f "$tmp_bam"
        else
            mv "$tmp_bam" {output.bam}
        fi
        samtools index {output.bam} {output.bai} >> {log} 2>&1
        """
