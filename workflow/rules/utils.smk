rule bam_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bai",
    log:
        "<logs>/bam-index/{prefix}.log",
    wrapper:
        "v2.3.2/bio/samtools/index"


rule tabix_known_variants:
    input:
        "<resources>/{prefix}.{format}.gz",
    output:
        "<resources>/{prefix}.{format}.gz.tbi",
    log:
        "<logs>/tabix/{prefix}.{format}.log",
    params:
        extra=branch(
            evaluate("{format}"),
            cases={
                "vcf": "-p vcf",
                "txt": "-s 1 -b 2 -e 2",
            },
        ),
    cache: "omit-software"
    wrapper:
        "v2.3.2/bio/tabix/index"


rule get_final_bam:
    input:
        bam=get_final_bam_input,
    output:
        bam="<results>/{sample}.bam",
    log:
        "<logs>/final-bam/{sample}.log",
    shell:
        "ln {input.bam} {output.bam} 2> {log}"
