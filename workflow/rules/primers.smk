rule assign_primers:
    input:
        bam=get_trimming_input,
        primers=get_primer_regions,
    output:
        assigned=temp("<results>/primers/{sample}.assigned.bam"),
        metric="<results>/primers/{sample}.metric.bam",
    log:
        "<logs>/primers/assignment/{sample}.log",
    conda:
        "../envs/fgbio.yaml"
    shell:
        "fgbio AssignPrimers -i {input.bam} -p {input.primers} -m {output.metric} -o {output.assigned} &> {log}"


rule filter_primerless_reads:
    input:
        "<results>/primers/{sample}.assigned.bam",
    output:
        primers=temp("<results>/primers/{sample}.primers.bam"),
        primerless=temp("<results>/primers/{sample}.primerless.bam"),
    log:
        "<logs>/primers/filter/{sample}.log",
    conda:
        "../envs/filter_reads.yaml"
    script:
        "../scripts/filter_primers.rs"


rule trim_primers:
    input:
        bam="<results>/primers/{sample}.primers.bam",
        primers=get_primer_regions,
    output:
        trimmed=branch(
            trimmed_is_final(),
            then="<results_trimmed>/{sample}.bam",
            otherwise=temp("<results_trimmed>/{sample}.bam"),
        ),
    log:
        "<logs>/trimming/{sample}.log",
    conda:
        "../envs/fgbio.yaml"
    params:
        sort_order="Coordinate",
        single_primer=branch(
            lambda wc: not isinstance(get_sample_primer_fastas(wc.sample), list),
            then="--first-of-pair",
            otherwise="",
        ),
    shell:
        "fgbio TrimPrimers -H -i {input.bam} -p {input.primers} -s {params.sort_order} {params.single_primer} -o {output.trimmed} &> {log}"


rule map_primers:
    input:
        reads=lambda wc: get_panel_primer_input(wc.panel),
        idx=access.random(rules.bwa_index.output),
    output:
        "<results>/primers/{panel}_primers.bam",
    log:
        "<logs>/bwa_mem/{panel}.log",
    threads: 8
    params:
        extra=get_primer_extra,
        sorting="none",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
    wrapper:
        "v2.13.0/bio/bwa/mem"


rule filter_unmapped_primers:
    input:
        "<results>/primers/{panel}_primers.bam",
    output:
        "<results>/primers/{panel}_primers.filtered.bam",
    log:
        "<logs>/primers/{panel}_primers_filtered.log",
    params:
        extra=branch(
            lambda wc: isinstance(get_panel_primer_input(wc.panel), list),
            then="-b -F 12",
            otherwise="-b -F 4",
        ),
    wrapper:
        "v2.3.2/bio/samtools/view"


rule primer_to_bed:
    input:
        "<results>/primers/{panel}_primers.filtered.bam",
    output:
        "<results>/primers/{panel}_primers.{ext}",
    log:
        "<logs>/primers/{panel}_primers_{ext}.log",
    wildcard_constraints:
        ext="bedpe|bed",
    conda:
        "../envs/bedtools.yaml"
    params:
        format=branch(
            lambda wc: wc.ext == "bedpe",
            then="-bedpe",
            otherwise="",
        ),
    shell:
        "samtools sort -n {input} | bamToBed -i - {params.format} > {output} 2> {log}"


rule build_primer_regions:
    input:
        branch(
            isinstance(primer_panels, pd.DataFrame),
            then=branch(
                lambda wc: not pd.isna(primer_panels.loc[wc.panel, "fa2"]),
                then="<results>/primers/{panel}_primers.bedpe",
                otherwise="<results>/primers/{panel}_primers.bed",
            ),
            otherwise=branch(
                bool(config["primers"]["trimming"].get("primers_fa2", "")),
                then="<results>/primers/uniform_primers.bedpe",
                otherwise="<results>/primers/uniform_primers.bed",
            ),
        ),
    output:
        "<results>/primers/{panel}_primer_regions.tsv",
    log:
        "<logs>/primers/build_{panel}_primer_regions.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_primer_regions.py"
