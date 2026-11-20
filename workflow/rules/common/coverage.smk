from clip import (
    coverage_normalization_mode,
    representative_count_samtools_args,
    track_strand_args,
)
from params import module_tool_extra, render_tool_params, resource_value
from tracks import coverage_strand_arg


if module_settings("coverage").get("strand_mode") == "separate":
    rule bam_coverage_stranded:
        input:
            bam=lambda wildcards: final_bam_path(wildcards.sample),
            bai=lambda wildcards: final_bai_path(wildcards.sample)
        output:
            bw=f"{RESULTS_DIR}/tracks/read_coverage/{{sample}}.{{strand}}.bw"
        wildcard_constraints:
            strand="plus|minus"
        log:
            f"{RESULTS_DIR}/logs/deeptools/read_coverage.{{sample}}.{{strand}}.log"
        threads:
            int(resource_value(config, "deeptools", "threads", 4))
        resources:
            mem_mb=int(resource_value(config, "deeptools", "mem_mb", 4096)),
            runtime_min=int(resource_value(config, "deeptools", "runtime_min", 120))
        conda:
            workflow_file("envs/deeptools.yaml")
        params:
            norm_mode=lambda wildcards: coverage_normalization_mode(config),
            count_args=lambda wildcards: representative_count_samtools_args(
                config, SAMPLES, wildcards.sample, "coverage"
            ),
            rendered=lambda wildcards: render_tool_params(
                config,
                "deeptools",
                module="coverage",
                overrides={"normalize_using": "None"},
                drop_keys=("effective_genome_size",),
            ),
            strand=lambda wildcards: track_strand_args(
                config, SAMPLES, wildcards.sample, wildcards.strand, "coverage"
            ),
            extra=lambda wildcards: module_tool_extra(config, "coverage", "deeptools")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.bw:q}) $(dirname {log:q})
            scale_factor="1"
            if [ "{params.norm_mode}" = "CPM" ]; then
                total_reads=$(samtools view -@ {threads} -c {params.count_args} {input.bam:q})
                scale_factor=$(python -c 'import sys; n = int(sys.argv[1]); print(f"{{1000000.0 / n:.12g}}" if n > 0 else "1")' "$total_reads")
            fi
            bamCoverage -b {input.bam:q} -o {output.bw:q} -p {threads} \
                --scaleFactor "$scale_factor" \
                {params.rendered} {params.strand} {params.extra} > {log:q} 2>&1
            """
else:
    rule bam_coverage:
        input:
            bam=lambda wildcards: final_bam_path(wildcards.sample),
            bai=lambda wildcards: final_bai_path(wildcards.sample)
        output:
            bw=f"{RESULTS_DIR}/tracks/bigwig/{{sample}}.bw"
        log:
            f"{RESULTS_DIR}/logs/deeptools/bamCoverage.{{sample}}.log"
        threads:
            int(resource_value(config, "deeptools", "threads", 4))
        resources:
            mem_mb=int(resource_value(config, "deeptools", "mem_mb", 4096)),
            runtime_min=int(resource_value(config, "deeptools", "runtime_min", 120))
        conda:
            workflow_file("envs/deeptools.yaml")
        params:
            rendered=lambda wildcards: render_tool_params(
                config,
                "deeptools",
                module="coverage",
                overrides={
                    "effective_genome_size": config.get("reference", {}).get(
                        "effective_genome_size"
                    )
                },
            ),
            strand=lambda wildcards: coverage_strand_arg(config, SAMPLES, wildcards.sample),
            extra=lambda wildcards: module_tool_extra(config, "coverage", "deeptools")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.bw}) $(dirname {log})
            bamCoverage -b {input.bam} -o {output.bw} -p {threads} \
                {params.rendered} {params.strand} {params.extra} > {log} 2>&1
            """
