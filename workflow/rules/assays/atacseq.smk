from params import render_tool_params, resource_value, tool_extra
from samples import (
    fragment_length_targets,
    frip_targets,
    library_complexity_targets,
    narrowpeak_targets,
    sample_layout,
    tss_enrichment_targets,
)


ATACSEQ_TARGETS = []

if output_enabled("fastq_qc"):
    ATACSEQ_TARGETS.extend(fastqc_targets(SAMPLES, RESULTS_DIR))

if output_enabled("trimmed_fastq") and step_enabled("trimming", False):
    ATACSEQ_TARGETS.extend(trimmed_fastq_targets(SAMPLES, RESULTS_DIR))

if output_enabled("raw_bam"):
    ATACSEQ_TARGETS.extend(raw_bam_targets(SAMPLES, RESULTS_DIR))

if output_enabled("filtered_bam") and step_enabled("bam_filter", True):
    ATACSEQ_TARGETS.extend(filtered_bam_targets(SAMPLES, RESULTS_DIR))

if output_enabled("bam_qc") and step_enabled("bam_qc", True):
    ATACSEQ_TARGETS.extend(bam_qc_targets(SAMPLES, RESULTS_DIR))

if output_enabled("bigwig") and step_enabled("coverage", False):
    ATACSEQ_TARGETS.extend(bigwig_targets(SAMPLES, RESULTS_DIR))

if output_enabled("peaks") and step_enabled("peak_calling", True):
    ATACSEQ_TARGETS.extend(narrowpeak_targets(SAMPLES, RESULTS_DIR))

if output_enabled("frip") and step_enabled("frip", True):
    ATACSEQ_TARGETS.extend(frip_targets(SAMPLES, RESULTS_DIR))

if output_enabled("library_complexity") and step_enabled("library_complexity", True):
    ATACSEQ_TARGETS.extend(library_complexity_targets(SAMPLES, RESULTS_DIR))

if output_enabled("fragment_lengths") and step_enabled("fragment_lengths", True):
    ATACSEQ_TARGETS.extend(fragment_length_targets(SAMPLES, RESULTS_DIR))

if output_enabled("tss_enrichment") and step_enabled("tss_enrichment", True):
    ATACSEQ_TARGETS.extend(tss_enrichment_targets(SAMPLES, RESULTS_DIR))

MULTIQC_INPUTS.extend(ATACSEQ_TARGETS)

if output_enabled("multiqc"):
    ATACSEQ_TARGETS.append(f"{RESULTS_DIR}/multiqc/multiqc_report.html")

ASSAY_TARGETS.extend(ATACSEQ_TARGETS)


rule macs3_callpeak:
    input:
        bam=lambda wildcards: final_bam_path(wildcards.sample),
        bai=lambda wildcards: final_bai_path(wildcards.sample)
    output:
        narrowpeak=f"{RESULTS_DIR}/peaks/macs3/{{sample}}/{{sample}}_peaks.narrowPeak",
        xls=f"{RESULTS_DIR}/peaks/macs3/{{sample}}/{{sample}}_peaks.xls",
        summits=f"{RESULTS_DIR}/peaks/macs3/{{sample}}/{{sample}}_summits.bed"
    log:
        f"{RESULTS_DIR}/logs/macs3/{{sample}}.log"
    threads:
        int(resource_value(config, "macs3", "threads", 1))
    resources:
        mem_mb=int(resource_value(config, "macs3", "mem_mb", 4096)),
        runtime_min=int(resource_value(config, "macs3", "runtime_min", 120))
    conda:
        str(WORKFLOW_DIR / "envs" / "macs3.yaml")
    params:
        outdir=lambda wildcards, output: Path(output.narrowpeak).parent.as_posix(),
        name=lambda wildcards: wildcards.sample,
        fmt=lambda wildcards: "BAMPE" if sample_layout(SAMPLES, wildcards.sample) == "paired" else "BAM",
        rendered=lambda wildcards: render_tool_params(
            config,
            "macs3",
            section="callpeak",
            overrides={
                "genome_size": config.get("tools", {}).get("macs3", {}).get("params", {}).get("callpeak", {}).get("genome_size")
                or config.get("genome", {}).get("macs3_genome_size")
                or config.get("genome", {}).get("effective_genome_size")
            },
        ),
        extra=lambda wildcards: tool_extra(config, "macs3")
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        macs3 callpeak -t {input.bam} -f {params.fmt} -n {params.name} \
            --outdir {params.outdir} {params.rendered} {params.extra} > {log} 2>&1
        """


rule atac_frip:
    input:
        bam=lambda wildcards: final_bam_path(wildcards.sample),
        bai=lambda wildcards: final_bai_path(wildcards.sample),
        peaks=f"{RESULTS_DIR}/peaks/macs3/{{sample}}/{{sample}}_peaks.narrowPeak"
    output:
        f"{RESULTS_DIR}/qc/atac/{{sample}}.frip.txt"
    log:
        f"{RESULTS_DIR}/logs/atac/frip.{{sample}}.log"
    threads:
        int(resource_value(config, "atac_frip", "threads", 1))
    resources:
        mem_mb=int(resource_value(config, "atac_frip", "mem_mb", 1024)),
        runtime_min=int(resource_value(config, "atac_frip", "runtime_min", 30))
    conda:
        str(WORKFLOW_DIR / "envs" / "bedtools.yaml")
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        total=$(samtools view -@ {threads} -c {input.bam})
        in_peaks=$(bedtools intersect -u -abam {input.bam} -b {input.peaks} \
            | samtools view -@ {threads} -c -)
        awk -v total="$total" -v in_peaks="$in_peaks" 'BEGIN {{
            frip = (total > 0) ? in_peaks / total : 0;
            print "sample_id\\ttotal_reads\\treads_in_peaks\\tfrip";
            print "{wildcards.sample}\\t" total "\\t" in_peaks "\\t" frip;
        }}' > {output} 2> {log}
        """


rule atac_library_complexity:
    input:
        bam=lambda wildcards: raw_bam_path(wildcards.sample),
        bai=lambda wildcards: raw_bai_path(wildcards.sample)
    output:
        f"{RESULTS_DIR}/qc/atac/{{sample}}.library_complexity.txt"
    log:
        f"{RESULTS_DIR}/logs/atac/library_complexity.{{sample}}.log"
    threads:
        int(resource_value(config, "atac_library_complexity", "threads", 1))
    resources:
        mem_mb=int(resource_value(config, "atac_library_complexity", "mem_mb", 1024)),
        runtime_min=int(resource_value(config, "atac_library_complexity", "runtime_min", 30))
    conda:
        str(WORKFLOW_DIR / "envs" / "samtools.yaml")
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        samtools view -@ {threads} {input.bam} \
            | awk 'BEGIN {{ OFS="\\t" }}
                {{
                    strand = (int($2 / 16) % 2) ? "-" : "+";
                    key = $3 ":" $4 ":" strand;
                    count[key]++;
                    total++;
                }}
                END {{
                    for (key in count) {{
                        distinct++;
                        if (count[key] == 1) one++;
                        if (count[key] == 2) two++;
                    }}
                    nrf = (total > 0) ? distinct / total : 0;
                    pbc1 = (distinct > 0) ? one / distinct : 0;
                    pbc2 = (two > 0) ? one / two : "Inf";
                    print "sample_id", "total_reads", "distinct_reads", "one_read_positions", "two_read_positions", "NRF", "PBC1", "PBC2";
                    print "{wildcards.sample}", total, distinct, one, two, nrf, pbc1, pbc2;
                }}' > {output} 2> {log}
        """


rule atac_fragment_lengths:
    input:
        bam=lambda wildcards: final_bam_path(wildcards.sample),
        bai=lambda wildcards: final_bai_path(wildcards.sample)
    output:
        f"{RESULTS_DIR}/qc/atac/{{sample}}.fragment_lengths.txt"
    log:
        f"{RESULTS_DIR}/logs/atac/fragment_lengths.{{sample}}.log"
    threads:
        int(resource_value(config, "atac_fragment_lengths", "threads", 1))
    resources:
        mem_mb=int(resource_value(config, "atac_fragment_lengths", "mem_mb", 1024)),
        runtime_min=int(resource_value(config, "atac_fragment_lengths", "runtime_min", 30))
    conda:
        str(WORKFLOW_DIR / "envs" / "samtools.yaml")
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        samtools stats -@ {threads} {input.bam} > {log} 2>&1
        awk 'BEGIN {{ OFS="\\t"; print "insert_size", "count" }} $1 == "IS" {{ print $2, $3 }}' {log} > {output}
        """


rule atac_tss_enrichment:
    input:
        bam=lambda wildcards: final_bam_path(wildcards.sample),
        bai=lambda wildcards: final_bai_path(wildcards.sample),
        tss=lambda wildcards: config.get("genome", {}).get("tss_bed"),
        chrom_sizes=lambda wildcards: config.get("genome", {}).get("chrom_sizes")
    output:
        f"{RESULTS_DIR}/qc/atac/{{sample}}.tss_enrichment.txt"
    log:
        f"{RESULTS_DIR}/logs/atac/tss_enrichment.{{sample}}.log"
    threads:
        int(resource_value(config, "atac_tss_enrichment", "threads", 1))
    resources:
        mem_mb=int(resource_value(config, "atac_tss_enrichment", "mem_mb", 1024)),
        runtime_min=int(resource_value(config, "atac_tss_enrichment", "runtime_min", 30))
    conda:
        str(WORKFLOW_DIR / "envs" / "bedtools.yaml")
    params:
        flank=lambda wildcards: int(config.get("atac", {}).get("tss_flank_bp", 1000)),
        center=lambda wildcards: int(config.get("atac", {}).get("tss_center_bp", 50))
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output}) $(dirname {log})
        tmpdir=$(mktemp -d)
        bedtools slop -i {input.tss} -g {input.chrom_sizes} -b {params.flank} > "$tmpdir/tss_window.bed"
        bedtools slop -i {input.tss} -g {input.chrom_sizes} -b {params.center} > "$tmpdir/tss_center.bed"
        window_reads=$(bedtools intersect -u -abam {input.bam} -b "$tmpdir/tss_window.bed" | samtools view -@ {threads} -c -)
        center_reads=$(bedtools intersect -u -abam {input.bam} -b "$tmpdir/tss_center.bed" | samtools view -@ {threads} -c -)
        awk -v window_reads="$window_reads" -v center_reads="$center_reads" -v flank="{params.flank}" -v center="{params.center}" 'BEGIN {{
            expected = (window_reads > 0 && flank > 0) ? window_reads * ((2 * center) / (2 * flank)) : 0;
            enrichment = (expected > 0) ? center_reads / expected : 0;
            print "sample_id\\twindow_reads\\tcenter_reads\\ttss_enrichment_proxy";
            print "{wildcards.sample}\\t" window_reads "\\t" center_reads "\\t" enrichment;
        }}' > {output} 2> {log}
        rm -rf "$tmpdir"
        """
