rule bam_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bai",
    log:
        "results/logs/bam-index/{prefix}.log",
    wrapper:
        "v8.1.1/bio/samtools/index"