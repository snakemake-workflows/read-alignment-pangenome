rule map_reads_bwa:
    input:
        reads=get_map_reads_input,
        idx=access.random(rules.bwa_index.output),
    output:
        temp("results/mapped/bwa/{sample}.bam"),
    log:
        "results/logs/bwa_mem/{sample}.log",
    params:
        extra=get_read_group("-R "),
    threads: 8
    wrapper:
        "v8.1.1/bio/bwa/mem"


rule count_sample_kmers:
    input:
        reads=get_map_reads_input,
    output:
        "results/kmers/{sample}.kff",
    params:
        out_file=lambda wc, output: os.path.splitext(output[0])[0],
        mem=lambda wc, resources: int(resources.mem_gb),
    conda:
        "../envs/kmc.yaml"
    shadow:
        "minimal"
    log:
        "results/logs/kmers/{sample}.log",
    threads: min(int(config.get("kmc_threads", 8)), 128)
    resources:
        mem_gb=64,
    shell:
        "kmc -k29 -m{params.mem} -sm -okff -t{threads} -v @<(printf '%s\n' {input.reads}) "
        '"{params.out_file}" . &> {log}'


rule create_reference_paths:
    output:
        "resources/reference_paths.txt",
    params:
        build=config["ref"]["build"],
    conda:
        "../envs/coreutils.yaml"
    log:
        "results/logs/reference/paths.log",
    shell:
        'for chrom in {{1..22}} X Y M; do echo "{params.build}#0#chr$chrom"; done > {output} 2> {log}'


rule map_reads_vg:
    input:
        reads=get_map_reads_input,
        graph=f"{pangenome_prefix}.gbz",
        kmers="results/kmers/{sample}.kff",
        hapl=f"{pangenome_prefix}.hapl",
        paths="resources/reference_paths.txt",
    output:
        bam=temp("results/mapped/vg/{sample}.raw.bam"),
        indexes=temp(
            multiext(
                f"{pangenome_prefix}.{{sample}}",
                ".gbz",
                ".dist",
                ".shortread.withzip.min",
                ".shortread.zipcodes",
            )
        ),
    log:
        "results/logs/mapped/vg/{sample}.log",
    benchmark:
        "results/benchmarks/vg_giraffe/{sample}.tsv"
    params:
        extra=lambda wc, input: f"--ref-paths {input.paths}",
        sorting="none",
    threads: 64
    wrapper:
        "v8.1.1/bio/vg/giraffe"


rule reheader_mapped_reads:
    input:
        "results/mapped/vg/{sample}.raw.bam",
    output:
        temp("results/mapped/vg/{sample}.reheadered.bam"),
    params:
        build=config["ref"]["build"],
    conda:
        "../envs/samtools.yaml"
    log:
        "results/logs/reheader/{sample}.log",
    shell:
        "(samtools view {input} -H | "
        "sed -E 's/(SN:{params.build}#0#chr)/SN:/; s/SN:M(\\t|$)/SN:MT\\1/' | "
        "samtools reheader - {input} > {output}) 2> {log}"


rule fix_mate:
    input:
        "results/mapped/vg/{sample}.reheadered.bam",
    output:
        temp("results/mapped/vg/{sample}.mate_fixed.bam"),
    log:
        "results/logs/samtools/fix_mate/{sample}.log",
    threads: 8
    params:
        extra="",
    wrapper:
        "v8.1.1/bio/samtools/fixmate"


# adding read groups is exclusive to vg mapped reads and
# necessary because base recalibration throws errors
# for not being able to find read group information
rule add_read_group:
    input:
        lambda wc: (
            "results/mapped/vg/{sample}.mate_fixed.bam"
            if sample_has_primers(wc)
            else "results/mapped/vg/{sample}.reheadered.bam"
        ),
    output:
        temp("results/mapped/vg/{sample}.bam"),
    log:
        "results/logs/samtools/add_rg/{sample}.log",
    params:
        read_group=get_read_group(""),
        compression_threads=lambda wildcards, threads: (
            f"-@{threads}" if threads > 1 else ""
        ),
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools addreplacerg {input} -o {output} -r {params.read_group} "
        "-w {params.compression_threads} 2> {log}"


rule sort_alignments:
    input:
        "results/mapped/{aligner}/{sample}.bam",
    output:
        temp("results/mapped/{aligner}/{sample}.sorted.bam"),
    log:
        "results/logs/sort/{aligner}/{sample}.log",
    threads: 16
    resources:
        mem_mb=32000,
    wrapper:
        "v8.1.1/bio/samtools/sort"


rule annotate_umis:
    input:
        bam="results/mapped/{aligner}/{sample}.sorted.bam",
        idx="results/mapped/{aligner}/{sample}.sorted.bai",
    output:
        temp("results/mapped/{aligner}/{sample}.annotated.bam"),
    conda:
        "../envs/umi_tools.yaml"
    log:
        "results/logs/annotate_bam/{aligner}/{sample}.log",
    shell:
        "umi_tools group -I {input.bam} --paired --umi-separator : --output-bam -S {output} &> {log}"


rule mark_duplicates:
    input:
        bams=get_markduplicates_input,
    output:
        bam=temp("results/dedup/{sample}.bam"),
        metrics="results/qc/dedup/{sample}.metrics.txt",
    log:
        "results/logs/picard/dedup/{sample}.log",
    params:
        extra=get_markduplicates_extra,
    resources:
        mem_mb=3000,
    wrapper:
        "v8.1.1/bio/picard/markduplicates"


rule calc_consensus_reads:
    input:
        get_consensus_input,
    output:
        consensus_r1=temp("results/consensus/fastq/{sample}.1.fq"),
        consensus_r2=temp("results/consensus/fastq/{sample}.2.fq"),
        consensus_se=temp("results/consensus/fastq/{sample}.se.fq"),
        skipped=temp("results/consensus/{sample}.skipped.bam"),
    log:
        "results/logs/consensus/{sample}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt collapse-reads-to-fragments bam {input} {output} &> {log}"


rule map_consensus_reads:
    input:
        reads=get_processed_consensus_input,
        idx=access.random(rules.bwa_index.output),
    output:
        temp("results/consensus/{sample}.consensus.{read_type}.mapped.bam"),
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=lambda w: f"-C {get_read_group('-R')(w)}",
        sorting="samtools",
        sort_order="coordinate",
    wildcard_constraints:
        read_type="pe|se",
    log:
        "results/logs/bwa_mem/{sample}.{read_type}.consensus.log",
    threads: 8
    wrapper:
        "v8.1.1/bio/bwa/mem"


rule merge_consensus_reads:
    input:
        "results/consensus/{sample}.skipped.bam",
        "results/consensus/{sample}.consensus.se.mapped.bam",
        "results/consensus/{sample}.consensus.pe.mapped.bam",
    output:
        temp("results/consensus/{sample}.merged.bam"),
    log:
        "results/logs/samtools_merge/{sample}.log",
    params:
        extra="",
    threads: 8
    wrapper:
        "v8.1.1/bio/samtools/merge"


rule sort_consensus_reads:
    input:
        "results/consensus/{sample}.merged.bam",
    output:
        temp("results/consensus/{sample}.bam"),
    log:
        "results/logs/samtools_sort/{sample}.log",
    threads: 16
    resources:
        mem_mb=64000,
    wrapper:
        "v8.1.1/bio/samtools/sort"


rule splitncigarreads:
    input:
        bam=lambda wc: (
            "results/dedup/{sample}.bam"
            if is_activated("remove_duplicates")
            else "results/mapped/star/{sample}.bam"
        ),
        ref=genome,
    output:
        "results/split/{sample}.bam",
    log:
        "results/logs/gatk/splitNCIGARreads/{sample}.log",
    params:
        extra="",
        java_opts="",
    resources:
        mem_mb=1024,
    wrapper:
        "v8.1.1/bio/gatk/splitncigarreads"


rule recalibrate_base_qualities:
    input:
        bam=get_recalibrate_quality_input,
        bai=lambda w: get_recalibrate_quality_input(w, bai=True),
        ref=genome,
        ref_dict=genome_dict,
        ref_fai=genome_fai,
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
    output:
        recal_table=temp("results/recal/{sample}.grp"),
    params:
        extra=config["params"]["gatk"]["BaseRecalibrator"],
        java_opts="",
    resources:
        mem_mb=1024,
    log:
        "results/logs/gatk/baserecalibrator/{sample}.log",
    threads: 8
    wrapper:
        "v8.1.1/bio/gatk/baserecalibratorspark"


ruleorder: apply_bqsr > bam_index


rule apply_bqsr:
    input:
        bam=get_recalibrate_quality_input,
        bai=lambda w: get_recalibrate_quality_input(w, bai=True),
        ref=genome,
        ref_dict=genome_dict,
        ref_fai=genome_fai,
        recal_table="results/recal/{sample}.grp",
    output:
        bam=protected("results/recal/{sample}.bam"),
        bai="results/recal/{sample}.bai",
    log:
        "results/logs/gatk/gatk_applybqsr/{sample}.log",
    params:
        extra=config["params"]["gatk"]["applyBQSR"],
        java_opts="",
    wrapper:
        "v8.1.1/bio/gatk/applybqsr"
