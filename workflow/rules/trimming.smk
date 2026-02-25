rule get_sra:
    output:
        "resources/sra/{accession}_1.fastq.gz",
        "resources/sra/{accession}_2.fastq.gz",
    log:
        "results/logs/get-sra/{accession}.log",
    wrapper:
        "v7.6.0/bio/sra-tools/fasterq-dump"


rule fastp_pipe:
    input:
        get_fastp_pipe_input,
    output:
        pipe("results/pipe/fastp/{sample}/{unit}.{fq}.{ext}"),
    log:
        "results/logs/pipe-fastqs/fastp/{sample}-{unit}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 1
    conda:
        "../envs/coreutils.yaml"
    shell:
        "cat {input} > {output} 2> {log}"


rule fastp_se:
    input:
        sample=lambda wc: get_fastp_input(wc),
    output:
        trimmed=temp("results/trimmed/{sample}/{unit}.single.fastq.gz"),
        html="results/trimmed/{sample}/{unit}.single.qc.html",
        json="results/trimmed/{sample}/{unit}.single.json",
    log:
        "results/logs/fastp/se/{sample}_{unit}.log",
    params:
        adapters=get_fastp_adapters,
        extra=get_fastp_extra,
    threads: 8
    wrapper:
        "v7.6.0/bio/fastp"


rule fastp_pe:
    input:
        sample=lambda wc: get_fastp_input(wc),
    output:
        trimmed=[
            temp("results/trimmed/{sample}/{unit}_R1.fastq.gz"),
            temp("results/trimmed/{sample}/{unit}_R2.fastq.gz"),
        ],
        html="results/trimmed/{sample}/{unit}.paired.qc.html",
        json="results/trimmed/{sample}/{unit}.paired.json",
    log:
        "results/logs/fastp/pe/{sample}_{unit}.log",
    params:
        adapters=get_fastp_adapters,
        extra=get_fastp_extra,
    threads: 8
    wrapper:
        "v7.6.0/bio/fastp"


rule merge_trimmed_fastqs:
    input:
        get_trimmed_fastqs,
    output:
        "results/merged/{sample}_{read}.fastq.gz",
    log:
        "results/logs/merge-fastqs/trimmed/{sample}_{read}.log",
    wildcard_constraints:
        read="single|R1|R2",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "cat {input} > {output} 2> {log}"
