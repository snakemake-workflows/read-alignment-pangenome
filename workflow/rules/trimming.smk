rule get_sra:
    output:
        "sra/{accession}_1.fastq.gz",
        "sra/{accession}_2.fastq.gz",
    params:
        extra=lookup(within=config, dpath="params/get_sra/extra", default=""),
    log:
        "<logs>/get-sra/{accession}.log",
    wrapper:
        "v7.6.0/bio/sra-tools/fasterq-dump"


rule fastp_pipe:
    input:
        # TODO: use get_raw_reads directly
        get_raw_reads,
    output:
        pipe("pipe/fastp/{sample}/{unit}.{fq}.{ext}"),
    log:
        "<logs>/pipe-fastqs/fastp/{sample}-{unit}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 0  # this does not need CPU
    shell:
        "cat {input} > {output} 2> {log}"


rule fastp_se:
    input:
        # TODO: I think this should just work as:
        # sample=get_fastp_input
        sample=get_fastp_input,
    output:
        trimmed=temp("<results>/trimmed/{sample}/{unit}.single.fastq.gz"),
        html="<results>/trimmed/{sample}/{unit}.single.qc.html",
        json="<results>/trimmed/{sample}/{unit}.single.json",
    log:
        "<logs>/fastp/se/{sample}_{unit}.log",
    params:
        adapters=get_fastp_adapters,
        extra=get_fastp_extra,
    threads: 8
    wrapper:
        "v6.2.0/bio/fastp"


rule fastp_pe:
    input:
        sample=get_fastp_input,
    output:
        trimmed=[
            temp("<results>/trimmed/{sample}/{unit}_R1.fastq.gz"),
            temp("<results>/trimmed/{sample}/{unit}_R2.fastq.gz"),
        ],
        html="<results>/trimmed/{sample}/{unit}.paired.qc.html",
        json="<results>/trimmed/{sample}/{unit}.paired.json",
    log:
        "<logs>/fastp/pe/{sample}_{unit}.log",
    params:
        adapters=get_fastp_adapters,
        extra=get_fastp_extra,
    threads: 8
    wrapper:
        "v6.2.0/bio/fastp"


rule merge_trimmed_fastqs:
    input:
        # TODO: try turning into a branch() function right here
        branch(
            lambda wc: units.loc[wc.sample, "adapters"].notna().all(),
            then=lambda wc: expand(
                "<results>/trimmed/{sample}/{unit}_{read}.fastq.gz",
                unit=units.loc[wc.sample, "unit_name"],
                sample=wc.sample,
                read=wc.read,
            ),
            otherwise=lambda wc: [
                read
                for unit in units.loc[wc.sample, "unit_name"]
                for read in get_raw_reads(
                    wc.sample,
                    unit,
                    "fq1" if wc.read in {"R1", "single"} else "fq2",
                )
            ],
        ),
    output:
        "<results>/merged/{sample}_{read}.fastq.gz",
    log:
        "<logs>/merge-fastqs/trimmed/{sample}_{read}.log",
    wildcard_constraints:
        read="single|R1|R2",
    shell:
        "cat {input} > {output} 2> {log}"
