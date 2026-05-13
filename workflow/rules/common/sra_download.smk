from params import resource_value, tool_extra
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
        str(WORKFLOW_DIR / "envs" / "sra_tools.yaml")
    params:
        accession=lambda wildcards: sample_by_unit(SAMPLES, wildcards.unit).get("sra_id", ""),
        extra=lambda wildcards: tool_extra(config, "sra_tools")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.r1}) $(dirname {log})
        tmpdir=$(mktemp -d)
        fasterq-dump {params.accession} --split-files --threads {threads} \
            --outdir "$tmpdir" {params.extra} > {log} 2>&1
        if [ -s "$tmpdir/{params.accession}_1.fastq" ]; then
            gzip -c "$tmpdir/{params.accession}_1.fastq" > {output.r1}
        elif [ -s "$tmpdir/{params.accession}.fastq" ]; then
            gzip -c "$tmpdir/{params.accession}.fastq" > {output.r1}
        else
            echo "No R1 FASTQ produced for {params.accession}" >> {log}
            rm -rf "$tmpdir"
            exit 1
        fi
        if [ -s "$tmpdir/{params.accession}_2.fastq" ]; then
            gzip -c "$tmpdir/{params.accession}_2.fastq" > {output.r2}
        else
            : > {output.r2}
        fi
        rm -rf "$tmpdir"
        """
