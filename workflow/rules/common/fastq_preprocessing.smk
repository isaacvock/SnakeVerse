from clip import barcode_adapter_args, demultiplexed_barcode_files, umi_extract_args
from params import module_tool_extra, resource_value
from samples import (
    optional_preprocessing_fastq_for_read,
    pool_fastq_for_read,
    pool_source_row,
    preprocessing_fastq_for_read,
    processed_unit_layout,
    row_is_paired,
    sample_by_unit,
    unit_layout,
)


if demultiplexing_active(config):
    rule cutadapt_inline_demultiplex:
        input:
            r1=lambda wildcards: pool_fastq_for_read(SAMPLES, wildcards.pool, "R1", RESULTS_DIR),
            r2=lambda wildcards: (
                [pool_fastq_for_read(SAMPLES, wildcards.pool, "R2", RESULTS_DIR)]
                if row_is_paired(pool_source_row(SAMPLES, wildcards.pool))
                else []
            )
        output:
            demux=temp(directory(f"{RESULTS_DIR}/fastq/demultiplexed/pools/{{pool}}")),
            json=f"{RESULTS_DIR}/qc/cutadapt/demultiplex.{{pool}}.json"
        log:
            f"{RESULTS_DIR}/logs/cutadapt/demultiplex.{{pool}}.log"
        threads:
            int(resource_value(config, "cutadapt_demultiplex", "threads", 4))
        resources:
            mem_mb=int(resource_value(config, "cutadapt_demultiplex", "mem_mb", 4096)),
            runtime_min=int(resource_value(config, "cutadapt_demultiplex", "runtime_min", 120))
        conda:
            workflow_file("envs/cutadapt.yaml")
        params:
            paired=lambda wildcards: row_is_paired(pool_source_row(SAMPLES, wildcards.pool)),
            adapters=lambda wildcards: barcode_adapter_args(config, SAMPLES, wildcards.pool),
            extra=lambda wildcards: module_tool_extra(
                config, "inline_barcode_demultiplexing", "cutadapt"
            )
        shell:
            r"""
            set -euo pipefail
            mkdir -p {output.demux:q} "$(dirname {output.json:q})" "$(dirname {log:q})"
            if [ "{params.paired}" = "True" ]; then
                cutadapt {params.adapters} {params.extra} --cores {threads} \
                    --json {output.json:q} \
                    -o '{output.demux}/{{name}}_R1.fastq.gz' \
                    -p '{output.demux}/{{name}}_R2.fastq.gz' \
                    {input.r1:q} {input.r2:q} > {log:q} 2>&1
            else
                cutadapt {params.adapters} {params.extra} --cores {threads} \
                    --json {output.json:q} \
                    -o '{output.demux}/{{name}}_R1.fastq.gz' \
                    {input.r1:q} > {log:q} 2>&1
            fi
            """


    rule merge_inline_barcodes:
        input:
            demux=lambda wildcards: f"{RESULTS_DIR}/fastq/demultiplexed/pools/{sample_by_unit(SAMPLES, wildcards.unit)['pool_id']}"
        output:
            r1=temp(f"{RESULTS_DIR}/fastq/demultiplexed/{{unit}}_R1.fastq.gz"),
            r2=temp(f"{RESULTS_DIR}/fastq/demultiplexed/{{unit}}_R2.fastq.gz")
        log:
            f"{RESULTS_DIR}/logs/cutadapt/merge.{{unit}}.log"
        params:
            r1=lambda wildcards: demultiplexed_barcode_files(
                RESULTS_DIR, SAMPLES, wildcards.unit, "R1"
            ),
            r2=lambda wildcards: demultiplexed_barcode_files(
                RESULTS_DIR, SAMPLES, wildcards.unit, "R2"
            ),
            paired=lambda wildcards: unit_layout(SAMPLES, wildcards.unit) == "paired"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.r1}) $(dirname {log})
            cat {params.r1} > {output.r1} 2> {log}
            if [ "{params.paired}" = "True" ]; then
                cat {params.r2} > {output.r2} 2>> {log}
            else
                gzip -c /dev/null > {output.r2}
            fi
            """


if umi_extraction_active(config):
    rule umi_tools_extract:
        input:
            r1=lambda wildcards: preprocessing_fastq_for_read(
                config, SAMPLES, wildcards.unit, "R1", RESULTS_DIR
            ),
            r2=lambda wildcards: optional_preprocessing_fastq_for_read(
                config, SAMPLES, wildcards.unit, "R2", RESULTS_DIR
            )
        output:
            r1=temp(f"{RESULTS_DIR}/fastq/umi_extracted/{{unit}}_R1.fastq.gz")
                if module_enabled("trimming")
                else f"{RESULTS_DIR}/fastq/umi_extracted/{{unit}}_R1.fastq.gz",
            r2=temp(f"{RESULTS_DIR}/fastq/umi_extracted/{{unit}}_R2.fastq.gz")
                if module_enabled("trimming")
                else f"{RESULTS_DIR}/fastq/umi_extracted/{{unit}}_R2.fastq.gz"
        log:
            f"{RESULTS_DIR}/logs/umi_tools/extract.{{unit}}.log"
        threads: 1
        resources:
            mem_mb=int(resource_value(config, "umi_tools_extract", "mem_mb", 4096)),
            runtime_min=int(resource_value(config, "umi_tools_extract", "runtime_min", 120))
        conda:
            workflow_file("envs/umi_tools.yaml")
        params:
            paired=lambda wildcards: unit_layout(SAMPLES, wildcards.unit) == "paired",
            rendered=lambda wildcards: umi_extract_args(
                config, unit_layout(SAMPLES, wildcards.unit) == "paired"
            ),
            extra=lambda wildcards: module_tool_extra(config, "umi_extraction", "umi_tools")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.r1}) $(dirname {log})
            if [ "{params.paired}" = "True" ]; then
                umi_tools extract {params.rendered} {params.extra} \
                    --stdin={input.r1} --stdout={output.r1} \
                    --read2-in={input.r2} --read2-out={output.r2} \
                    --log={log}
            else
                umi_tools extract {params.rendered} {params.extra} \
                    --stdin={input.r1} --stdout={output.r1} --log={log}
                gzip -c /dev/null > {output.r2}
            fi
            """
