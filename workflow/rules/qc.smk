rule fastqc:
    input:
        get_fastqc_input,
    output:
        html="<results>/qc/fastqc/{sample}/{unit}.{fq}.html",
        zip="<results>/qc/fastqc/{sample}/{unit}.{fq}_fastqc.zip",
    log:
        "<results>/logs/fastqc/{sample}/{unit}.{fq}.log",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.10.0/bio/fastqc"

rule samtools_idxstats_cram:
    input:
        cram="<results>/mapped/vg/{sample}.sorted.cram",
        crai="<results>/mapped/vg/{sample}.sorted.cram.crai",
        ref=genome,
    output:
        "<results>/qc/{sample}.cram.idxstats",
    log:
        "<results>/logs/samtools/idxstats/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 1
    shell:
        "samtools idxstats {input.cram} > {output} 2> {log}"

rule samtools_stats_cram:
    input:
        cram="<results>/mapped/vg/{sample}.sorted.cram",
        ref=genome,
    output:
        "<results>/qc/{sample}.cram.stats",
    log:
        "<results>/logs/samtools/stats/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 2
    shell:
        "samtools stats -@ {threads} -T {input.ref} {input.cram} > {output} 2> {log}"

rule multiqc:
    input:
        get_fastqc_results,
    output:
        report(
            "<results>/qc/multiqc/{group}.html",
            category="Quality control",
            caption="../report/multiqc.rst",
            labels={"Sample group": "{group}"},
        ),
    params:
        "--exclude snippy",
    log:
        "<results>/logs/multiqc/{group}.log",
    wrapper:
        "v2.10.0/bio/multiqc"

