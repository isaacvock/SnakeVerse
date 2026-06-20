from params import render_tool_params, resource_value, tool_extra
from samples import sample_by_unit


rule sra_fastq:
    output:
        r1=f"{RESULTS_DIR}/fastq/sra/{{unit}}_R1.fastq.gz",
        r2=f"{RESULTS_DIR}/fastq/sra/{{unit}}_R2.fastq.gz"
    log:
        f"{RESULTS_DIR}/logs/sra/{{unit}}.log"
    threads:
        int(resource_value(config, "sra_tools", "threads", 4))
    resources:
        mem_mb=int(resource_value(config, "sra_tools", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "sra_tools", "runtime_min", 240))
    conda:
        workflow_file("envs/sra_tools.yaml")
    params:
        accession=lambda wildcards: sample_by_unit(SAMPLES, wildcards.unit).get("sra_id", ""),
        rendered=lambda wildcards: render_tool_params(config, "sra_tools", section="fasterq_dump"),
        extra=lambda wildcards: tool_extra(config, "sra_tools", section="fasterq_dump")
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.r1:q})" "$(dirname {log:q})"
        exec > {log:q} 2>&1

        tmpdir=$(mktemp -d)
        trap 'rm -rf "$tmpdir"' EXIT

        accession={params.accession:q}
        paired_r1="$tmpdir/$accession"_1.fastq
        paired_r2="$tmpdir/$accession"_2.fastq
        single_r1="$tmpdir/$accession".fastq

        fasterq-dump "$accession" --split-files --threads {threads} \
            --mem {resources.mem_mb}M --temp "$tmpdir" --outdir "$tmpdir" \
            {params.rendered} {params.extra}

        if [ -s "$paired_r1" ]; then
            pigz -p {threads} -c "$paired_r1" > {output.r1:q}
        elif [ -s "$single_r1" ]; then
            pigz -p {threads} -c "$single_r1" > {output.r1:q}
        else
            echo "No R1 FASTQ produced for $accession"
            exit 1
        fi

        if [ -s "$paired_r2" ]; then
            pigz -p {threads} -c "$paired_r2" > {output.r2:q}
        else
            pigz -p {threads} -c /dev/null > {output.r2:q}
        fi
        """
