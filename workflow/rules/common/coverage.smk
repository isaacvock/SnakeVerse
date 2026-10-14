from clip import track_strand_args
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
            strand=lambda wildcards: track_strand_args(
                config, SAMPLES, wildcards.sample, wildcards.strand, "coverage"
            ),
            extra=lambda wildcards: module_tool_extra(config, "coverage", "deeptools")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.bw}) $(dirname {log})
            if ! bamCoverage -b {input.bam} -o {output.bw} -p {threads} \
                {params.rendered} {params.strand} {params.extra} > {log} 2>&1; then
                if grep -q "ZeroDivisionError" {log}; then
                    echo "SnakeVerse: selected strand is empty; writing an unnormalized zero track." >> {log}
                    bamCoverage -b {input.bam} -o {output.bw} -p {threads} \
                        {params.rendered} {params.strand} {params.extra} \
                        --normalizeUsing None >> {log} 2>&1
                else
                    exit 1
                fi
            fi
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
            mkdir -p $(dirname {output.bw}) $(dirname {log})
            bamCoverage -b {input.bam} -o {output.bw} -p {threads} \
                {params.rendered} {params.strand} {params.extra} > {log} 2>&1
            """
