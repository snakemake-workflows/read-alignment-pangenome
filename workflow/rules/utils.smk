rule bam_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bai",
    log:
        "logs/bam-index/{prefix}.log",
    wrapper:
        "v2.3.2/bio/samtools/index"


rule tabix_known_variants:
    input:
        "<resources>/{prefix}.{format}.gz",
    output:
        "<resources>/{prefix}.{format}.gz.tbi",
    log:
        "logs/tabix/{prefix}.{format}.log",
    params:
        # TODO: turn this into a branch() function right here,
        # using the cases= entry
        # note: not sure whether the `otherwise=` will work with the `cases=`
        get_tabix_params,
    cache: "omit-software"
    wrapper:
        "v2.3.2/bio/tabix/index"
