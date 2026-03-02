rule tabix_known_variants:
    input:
        "resources/{prefix}.{format}.gz",
    output:
        "resources/{prefix}.{format}.gz.tbi",
    log:
        "results/logs/tabix/{prefix}.{format}.log",
    params:
        get_tabix_params,
    cache: "omit-software"
    wrapper:
        "v8.1.1/bio/tabix/index"


rule bam_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bai",
    log:
        "results/logs/bam-index/{prefix}.log",
    wrapper:
        "v8.1.1/bio/samtools/index"
