rule count_sample_kmers:
    input:
        reads=get_map_reads_input,
    output:
        "<results>/kmers/{sample}.kff",
    params:
        out_file=lambda wc, output: os.path.splitext(output[0])[0],
        mem=lambda wc, resources: int(resources.mem_gb),
    conda:
        "../envs/kmc.yaml"
    shadow:
        "minimal"
    log:
        "<results>/logs/kmers/{sample}.log",
    threads: min(int(config.get("kmc_threads", 8)), 128)
    resources:
        mem_gb=64,
    shell:
        "tmpdir=$(mktemp -d); "
        "trap 'rm -rf \"$tmpdir\"' EXIT; "
        "kmc -k29 -m{params.mem} -sm -okff -t{threads} -v {input.reads} "
        "\"{params.out_file}\" \"$tmpdir\" &> {log}"

rule create_reference_paths:
    output:
        "<resources>/reference_paths.txt",
    params:
        build=config["ref"]["build"],
    conda:
        "../envs/coreutils.yaml"
    log:
        "<results>/logs/reference/paths.log",
    shell:
        'for chrom in {{1..22}} X Y M; do echo "{params.build}#0#chr$chrom"; done > {output} 2> {log}'

rule map_reads_vg:
    input:
        reads=get_map_reads_input,
        graph=f"{pangenome_prefix}.gbz",
        kmers="<results>/kmers/{sample}.kff",
        hapl=f"{pangenome_prefix}.hapl",
        paths="<resources>/reference_paths.txt",
    output:
        bam=temp("<results>/mapped/vg/{sample}.raw.bam"),
    log:
        "<results>/logs/mapped/vg/{sample}.log",
    benchmark:
        "<results>/benchmarks/vg_giraffe/{sample}.tsv",
    params:
        extra=lambda wc, input: f"--ref-paths {input.paths}",
        sorting="none",
    threads: 64
    wrapper:
        "v6.1.0/bio/vg/giraffe"

rule reheader_mapped_reads:
    input:
        "<results>/mapped/vg/{sample}.raw.bam",
    output:
        temp("<results>/mapped/vg/{sample}.reheadered.bam"),
    params:
        build=config["ref"]["build"],
    conda:
        "../envs/samtools.yaml"
    log:
        "<results>/logs/reheader/{sample}.log",
    shell:
        "(samtools view {input} -H | "
        "sed -E 's/(SN:{params.build}#0#chr)/SN:/; s/SN:M/SN:MT/' | "
        "samtools reheader - {input} > {output}) 2> {log}"

rule fix_mate:
    input:
        "<results>/mapped/vg/{sample}.reheadered.bam",
    output:
        temp("<results>/mapped/vg/{sample}.mate_fixed.bam"),
    log:
        "<results>/logs/samtools/fix_mate/{sample}.log",
    threads: 8
    wrapper:
        "v8.1.1/bio/samtools/fixmate"

# adding read groups is exclusive to vg mapped reads and
# necessary because base recalibration throws errors
# for not being able to find read group information
rule add_read_group:
    input:
        "<results>/mapped/vg/{sample}.mate_fixed.bam",
    output:
        temp("<results>/mapped/vg/{sample}.bam"),
    log:
        "<results>/logs/samtools/add_rg/{sample}.log",
    params:
        read_group=get_read_group(""),
        compression_threads=lambda wildcards, threads: (
            f"-@ {threads}" if threads > 1 else ""
        ),
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools addreplacerg {params.compression_threads} {input} -o {output} -r {params.read_group} "
        "2> {log}"

rule sort_alignments:
    input:
        "<results>/mapped/vg/{sample}.bam",
    output:
        temp("<results>/mapped/vg/{sample}.sorted.bam"),
    log:
        "<results>/logs/sort/vg/{sample}.log",
    threads: 16
    resources:
        mem_mb=32000,
    wrapper:
        "v8.1.1/bio/samtools/sort"

rule bam_to_cram:
    input:
        bam="<results>/mapped/vg/{sample}.sorted.bam",
        ref=genome,
        fai=genome_fai,
    output:
        protected("<results>/mapped/vg/{sample}.sorted.cram"),
    log:
        "<results>/logs/samtools/bam_to_cram/vg/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools view -@ {threads} -C -T {input.ref} -o {output} {input.bam} 2> {log}"

rule cram_index:
    input:
        cram="<results>/mapped/vg/{sample}.sorted.cram",
        ref=genome,
        fai=genome_fai,
    output:
        "<results>/mapped/vg/{sample}.sorted.cram.crai",
    log:
        "<results>/logs/samtools/cram_index/vg/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 1
    shell:
        "samtools index -c {input.cram} 2> {log}"
