rule bam_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bai",
    log:
        "<results>/logs/bam-index/{prefix}.log",
    wrapper:
        "v2.3.2/bio/samtools/index"


