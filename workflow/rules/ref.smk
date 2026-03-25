rule get_genome:
    output:
        genome,
    log:
        "logs/get-genome.log",
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
        "logs/genome-faidx.log",
    cache: "omit-software"
    wrapper:
        "v2.3.2/bio/samtools/faidx"


rule genome_dict:
    input:
        genome,
    output:
        genome_dict,
    log:
        "logs/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: "omit-software"
    shell:
        "samtools dict {input} > {output} 2> {log} "


rule get_known_variants:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai=genome_fai,
    output:
        vcf="resources/variation.vcf.gz",
    log:
        "logs/get-known-variants.log",
    params:
        species=config["ref"]["species"],
        release=config["ref"]["release"],
        build=config["ref"]["build"],
        type="all",
        chromosome=config["ref"].get("chromosome"),
    cache: "omit-software"
    wrapper:
        "v7.5.0/bio/reference/ensembl-variation"


rule remove_iupac_codes:
    input:
        "resources/variation.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz",
    log:
        "logs/fix-iupac-alleles.log",
    conda:
        "../envs/rbt.yaml"
    cache: "omit-software"
    shell:
        "(rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}) 2> {log}"


rule bwa_index:
    input:
        genome,
    output:
        idx=multiext(genome, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    cache: True
    wrapper:
        "v2.3.2/bio/bwa/index"


rule get_pangenome_vcf:
    output:
        f"{pangenome_prefix}.wave.vcf.gz",
    params:
        url=get_pangenome_url,
    log:
        "<logs>/pangenome/vcf.log",
    cache: "omit-software"
    shell:
        "curl -L {params.url} -o {output} 2> {log}"


rule vg_autoindex_giraffe:
    input:
        ref=genome,
        vcf=f"{pangenome_prefix}.wave.vcf.gz",
    output:
        multiext(
            pangenome_prefix,
            ".dist",
            ".shortread.zipcodes",
            ".shortread.withzip.min",
            ".giraffe.gbz",
        ),
    log:
        "<logs>/pangenome/autoindex.log",
    params:
        extra="",
    threads: 8
    wrapper:
        "v8.0.3/bio/vg/autoindex"
