rule get_genome:
    output:
        genome,
    log:
        "<logs>/get-genome.log",
    cache: "omit-software"
    params:
        species=lookup(within=config, dpath="ref/species"),
        datatype="dna",
        build=lookup(within=config, dpath="ref/build"),
        release=lookup(within=config, dpath="ref/release"),
        chromosome=lookup(within=config, dpath="ref/chromosome", default=None),
    wrapper:
        "v7.3.0/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        genome,
    output:
        genome_fai,
    log:
        "<logs>/genome-faidx.log",
    cache: "omit-software"
    wrapper:
        "v2.3.2/bio/samtools/faidx"


rule genome_dict:
    input:
        genome,
    output:
        genome_dict,
    log:
        "<logs>/samtools/create_dict.log",
    cache: "omit-software"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools dict {input} > {output} 2> {log} "


rule get_known_variants:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai=genome_fai,
    output:
        vcf="<resources>/variation.vcf.gz",
    log:
        "<logs>/get-known-variants.log",
    cache: "omit-software"
    params:
        species=lookup(within=config, dpath="ref/species"),
        release=lookup(within=config, dpath="ref/release"),
        build=lookup(within=config, dpath="ref/build"),
        type="all",
        chromosome=lookup(within=config, dpath="ref/chromosome", default=None),
    wrapper:
        "v7.5.0/bio/reference/ensembl-variation"


rule remove_iupac_codes:
    input:
        "<resources>/variation.vcf.gz",
    output:
        "<resources>/variation.noiupac.vcf.gz",
    log:
        "<logs>/fix-iupac-alleles.log",
    cache: "omit-software"
    conda:
        "../envs/rbt.yaml"
    shell:
        "(rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}) 2> {log}"


rule bwa_index:
    input:
        genome,
    output:
        idx=multiext(genome, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "<logs>/bwa_index.log",
    cache: True
    wrapper:
        "v2.3.2/bio/bwa/index"


rule get_pangenome:
    output:
        f"{pangenome_prefix}.{{ext}}",
    log:
        "<logs>/pangenome/{ext}.log",
    wildcard_constraints:
        ext="hapl|gbz",
    cache: "omit-software"
    params:
        url=lambda wc: get_pangenome_url(wc.ext),
    shell:
        "curl -o {output} {params.url} 2> {log}"
