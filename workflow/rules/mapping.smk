rule map_reads_bwa:
    input:
        reads=get_map_reads_input,
        idx=access.random(rules.bwa_index.output),
    output:
        temp("<results>/mapped/bwa/{sample}.raw.bam"),
    log:
        "<logs>/bwa_mem/{sample}.log",
    threads: 8
    params:
        extra=get_read_group("-R "),
    wrapper:
        "v3.8.0/bio/bwa/mem"


rule count_sample_kmers:
    input:
        reads=get_map_reads_input,
    output:
        "<results>/kmers/{sample}.kff",
    log:
        "<logs>/kmers/{sample}.log",
    shadow:
        "minimal"
    conda:
        "../envs/kmc.yaml"
    threads: min(max(workflow.cores, 1), 128)  # kmc can use 128 threads at most
    resources:
        mem="64GB",
    params:
        out_file=lambda wc, output: os.path.splitext(output[0])[0],
        mem=lambda wc, resources: resources.mem[:-2],
    shell:
        "tmpdir=$(mktemp -d); "
        "kmc -k29 -m{params.mem} -sm -okff -t{threads} -v @<(ls {input.reads}) {params.out_file} $tmpdir &> {log} && "
        "rm -r $tmpdir || (rm -r $tmpdir && exit 1)"


rule create_reference_paths:
    output:
        "<resources>/reference_paths.txt",
    log:
        "<logs>/reference/paths.log",
    params:
        build=lookup(within=config, dpath="ref/build"),
    shell:
        'for chrom in {{1..22}} X Y M; do echo "{params.build}#0#chr$chrom"; done > {output} 2> {log}'


rule map_reads_vg:
    input:
        reads=get_map_reads_input,
        graph=access.random(f"{pangenome_prefix}.gbz"),
        kmers=access.random("<results>/kmers/{sample}.kff"),
        hapl=access.random(f"{pangenome_prefix}.hapl"),
        paths=access.random("<resources>/reference_paths.txt"),
    output:
        bam=temp("<results>/mapped/vg/{sample}.raw.bam"),
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
        "<logs>/mapped/vg/{sample}.log",
    benchmark:
        "<benchmarks>/vg_giraffe/{sample}.tsv"
    threads: 64
    params:
        extra=lambda wc, input: f"--ref-paths {input.paths}",
        sorting="none",
    wrapper:
        "v6.1.0/bio/vg/giraffe"


rule reheader_mapped_reads:
    input:
        "<results>/mapped/vg/{sample}.raw.bam",
    output:
        temp("<results>/mapped/vg/{sample}.reheadered.bam"),
    log:
        "<logs>/reheader/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        build=lookup(within=config, dpath="ref/build"),
    shell:
        "(samtools view {input} -H |"
        " sed -E 's/(SN:{params.build}#0#chr)/SN:/; s/SN:M/SN:MT/' | "
        " samtools reheader - {input} > {output}) 2> {log}"


rule fix_mate:
    input:
        "<results>/mapped/vg/{sample}.reheadered.bam",
    output:
        temp("<results>/mapped/vg/{sample}.mate_fixed.bam"),
    log:
        "<logs>/samtools/fix_mate/{sample}.log",
    threads: 8
    params:
        extra="",
    wrapper:
        "v4.7.2/bio/samtools/fixmate"


# adding read groups is exclusive to vg mapped reads and
# necessary because base recalibration throws errors
# for not being able to find read group information
rule add_read_group:
    input:
        branch(
            sample_has_primers,
            then="<results>/mapped/vg/{sample}.mate_fixed.bam",
            otherwise="<results>/mapped/vg/{sample}.reheadered.bam",
        ),
    output:
        temp("<results_mapped_vg>/{sample}.bam"),
    log:
        "<logs>/samtools/add_rg/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 4
    params:
        read_group=get_read_group(""),
    shell:
        "samtools addreplacerg {input} -o {output} -r {params.read_group} "
        "-w -@{threads} 2> {log}"


rule sort_alignments:
    input:
        "<results>/mapped/bwa/{sample}.raw.bam",
    output:
        temp("<results_mapped_bwa>/{sample}.bam"),
    log:
        "<logs>/sort/bwa/{sample}.log",
    threads: 16
    resources:
        mem_mb=32000,
    wrapper:
        "v8.1.1/bio/samtools/sort"


rule annotate_umis:
    input:
        bam=get_mapped_stage_input,
        idx=lambda wc: get_mapped_stage_input(wc, bai=True),
    output:
        temp(
            "<results>/mapped/{aligner}/{{sample}}.annotated.bam".format(
                aligner=get_aligner
            )
        ),
    log:
        "<logs>/annotate_bam/{aligner}/{{sample}}.log".format(aligner=get_aligner),
    conda:
        "../envs/umi_tools.yaml"
    shell:
        "umi_tools group -I {input.bam} --paired --umi-separator : --output-bam -S {output} &> {log}"


rule mark_duplicates:
    input:
        bams=branch(
            sample_has_umis,
            then="<results>/mapped/{aligner}/{{sample}}.annotated.bam".format(
                aligner=get_aligner
            ),
            otherwise=get_mapped_stage_input,
        ),
    output:
        bam=temp("<results_dedup>/{sample}.bam"),
        metrics="<results>/qc/dedup/{sample}.metrics.txt",
    log:
        "<logs>/picard/dedup/{sample}.log",
    resources:
        mem_mb=3000,
    params:
        extra=get_markduplicates_extra,
    wrapper:
        "v2.5.0/bio/picard/markduplicates"


rule calc_consensus_reads:
    input:
        get_consensus_input,
    output:
        consensus_r1=temp("<results>/consensus/fastq/{sample}.1.fq"),
        consensus_r2=temp("<results>/consensus/fastq/{sample}.2.fq"),
        consensus_se=temp("<results>/consensus/fastq/{sample}.se.fq"),
        skipped=temp("<results>/consensus/{sample}.skipped.bam"),
    log:
        "<logs>/consensus/{sample}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt collapse-reads-to-fragments bam {input} {output} &> {log}"


rule map_consensus_reads:
    input:
        reads=branch(
            evaluate("{read_type} == 'se'"),
            then="<results>/consensus/fastq/{sample}.se.fq",
            otherwise=[
                "<results>/consensus/fastq/{sample}.1.fq",
                "<results>/consensus/fastq/{sample}.2.fq",
            ],
        ),
        idx=access.random(rules.bwa_index.output),
    output:
        temp("<results>/consensus/{sample}.consensus.{read_type}.mapped.bam"),
    log:
        "<logs>/bwa_mem/{sample}.{read_type}.consensus.log",
    wildcard_constraints:
        read_type="pe|se",
    threads: 8
    params:
        index=lambda wc, input: os.path.splitext(input.idx[0])[0],
        extra=lambda wc: f"-C {get_read_group('-R')(wc)}",
        sort="samtools",
        sort_order="coordinate",
    wrapper:
        "v2.3.2/bio/bwa/mem"


rule merge_consensus_reads:
    input:
        "<results>/consensus/{sample}.skipped.bam",
        "<results>/consensus/{sample}.consensus.se.mapped.bam",
        "<results>/consensus/{sample}.consensus.pe.mapped.bam",
    output:
        temp("<results>/consensus/{sample}.merged.bam"),
    log:
        "<logs>/samtools_merge/{sample}.log",
    threads: 8
    wrapper:
        "v2.3.2/bio/samtools/merge"


rule sort_consensus_reads:
    input:
        "<results>/consensus/{sample}.merged.bam",
    output:
        temp("<results_consensus>/{sample}.bam"),
    log:
        "<logs>/samtools_sort/{sample}.log",
    threads: 16
    resources:
        mem_mb=64000,
    wrapper:
        "v8.1.1/bio/samtools/sort"


rule recalibrate_base_qualities:
    input:
        bam=get_recalibrate_quality_input,
        bai=lambda wc: get_recalibrate_quality_input(wc, bai=True),
        ref=genome,
        ref_dict=genome_dict,
        ref_fai=genome_fai,
        known="<resources>/variation.noiupac.vcf.gz",
        tbi="<resources>/variation.noiupac.vcf.gz.tbi",
    output:
        recal_table=temp("<results>/recal/{sample}.grp"),
    log:
        "<logs>/gatk/baserecalibrator/{sample}.log",
    threads: 8
    resources:
        mem_mb=1024,
    params:
        extra=lookup(within=config, dpath="params/gatk/BaseRecalibrator"),
        java_opts="",
    wrapper:
        "v1.25.0/bio/gatk/baserecalibratorspark"


ruleorder: apply_bqsr > bam_index


rule apply_bqsr:
    input:
        bam=get_recalibrate_quality_input,
        bai=lambda wc: get_recalibrate_quality_input(wc, bai=True),
        ref=genome,
        ref_dict=genome_dict,
        ref_fai=genome_fai,
        recal_table="<results>/recal/{sample}.grp",
    output:
        bam=protected("<results_bqsr>/{sample}.bam"),
        bai="<results_bqsr>/{sample}.bai",
    log:
        "<logs>/gatk/gatk_applybqsr/{sample}.log",
    params:
        extra=lookup(within=config, dpath="params/gatk/applyBQSR"),
        java_opts="",
    wrapper:
        "v2.3.2/bio/gatk/applybqsr"
