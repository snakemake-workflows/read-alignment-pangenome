rule get_genome:
    output:
        genome,
    log:
        "results/logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        chromosome=config["ref"].get("chromosome"),
    cache: "omit-software"
    wrapper:
        "v7.3.0/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        genome,
    output:
        genome_fai,
    log:
        "results/logs/genome-faidx.log",
    cache: "omit-software"
    wrapper:
        "v8.1.1/bio/samtools/faidx"


rule genome_dict:
    input:
        genome,
    output:
        genome_dict,
    log:
        "results/logs/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: "omit-software"
    shell:
        "samtools dict {input} > {output} 2> {log}"


rule bwa_index:
    input:
        genome,
    output:
        idx=multiext(genome, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "results/logs/bwa_index.log",
    cache: True
    wrapper:
        "v8.1.1/bio/bwa/index"


rule get_pangenome:
    output:
        f"{pangenome_prefix}.{{ext}}",
    params:
        url=lambda wc: get_pangenome_url(wc.ext),
    wildcard_constraints:
        ext="hapl|gbz",
    log:
        "results/logs/pangenome/{ext}.log",
    cache: "omit-software"
    conda:
        "../envs/curl.yaml"
    shell:
        "curl -o {output} {params.url} 2> {log}"
