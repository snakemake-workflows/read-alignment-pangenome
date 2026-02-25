import glob
import os
from os import path

import pandas as pd
from snakemake.utils import validate
from snakemake.exceptions import WorkflowError

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype={"sample_name": str, "group": str},
        comment="#",
    )
    .set_index("sample_name", drop=False)
    .sort_index()
)

# construct genome name
datatype_genome = "dna"
species = config["ref"]["species"]
build = config["ref"]["build"]
release = config["ref"]["release"]
genome_name = f"genome.{datatype_genome}.{species}.{build}.{release}"
genome_prefix = f"resources/{genome_name}"
genome = f"{genome_prefix}.fasta"
genome_fai = f"{genome}.fai"
genome_dict = f"{genome_prefix}.dict"
pangenome_name = f"pangenome.{species}.{build}"
pangenome_prefix = f"resources/{pangenome_name}"

# cram variables (mini-workflow: CRAM is mandatory)
use_cram = True


def _group_or_sample(row):
    group = row.get("group", None)
    if pd.isnull(group):
        return row["sample_name"]
    return group


samples["group"] = [_group_or_sample(row) for _, row in samples.iterrows()]


def get_group_samples(group):
    return samples.loc[samples["group"] == group]["sample_name"]


units = (
    pd.read_csv(
        config["units"],
        sep="\t",
        dtype={"sample_name": str, "unit_name": str},
        comment="#",
    )
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)

validate(units, schema="../schemas/units.schema.yaml")

primer_panels = (
    pd.read_csv(
        config["primers"]["trimming"]["tsv"],
        sep="\t",
        dtype={"panel": str, "fa1": str, "fa2": str},
        comment="#",
    )
    .set_index(["panel"], drop=False)
    .sort_index()
)

if primer_panels.empty:
    raise WorkflowError(
        "Primers TSV is empty: config['primers']['trimming']['tsv']="
        f"{config['primers']['trimming']['tsv']!r}. "
        "Downstream helpers expect primer entries; an empty TSV will later cause "
        "KeyError when dereferencing primers_fa1. Please provide a non-empty primers TSV."
    )


def get_fastp_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        return get_sra_reads(wildcards.sample, wildcards.unit, ["1", "2"])

    fq1 = get_raw_reads(unit.sample_name, unit.unit_name, "fq1")

    if len(fq1) == 1:

        def get_reads(fq):
            return get_raw_reads(unit.sample_name, unit.unit_name, fq)[0]

    else:
        ending = ".gz" if unit["fq1"].endswith("gz") else ""

        def get_reads(fq):
            return f"results/pipe/fastp/{unit.sample_name}/{unit.unit_name}.{fq}.fastq{ending}"

    if pd.isna(unit["fq2"]):
        # single end sample
        return get_reads("fq1")
    else:
        # paired end sample
        return [get_reads("fq1"), get_reads("fq2")]


def get_sra_reads(sample, unit, fq):
    unit = units.loc[sample].loc[unit]
    # SRA sample (always paired-end for now)
    accession = unit["sra"]
    return expand(
        "resources/sra/{accession}_{read}.fastq.gz", accession=accession, read=fq
    )


def get_raw_reads(sample, unit, fq):
    pattern = units.loc[sample].loc[unit, fq]

    if pd.isna(pattern):
        assert fq.startswith("fq")
        fq = fq[len("fq") :]
        return get_sra_reads(sample, unit, fq)

    if type(pattern) is not str and len(pattern) > 1:
        raise ValueError(
            f"Multiple units.tsv entries found for sample '{sample}' and "
            f"unit '{unit}'.\n"
            "The units.tsv should contain only one entry for each combination "
            "of sample and unit.\n"
            "Found:\n"
            f"{pattern}"
        )

    if "*" in pattern:
        files = sorted(glob.glob(units.loc[sample].loc[unit, fq]))
        if not files:
            raise ValueError(
                "No raw fastq files found for unit pattern {} (sample {}). "
                "Please check your sample sheet.".format(unit, sample)
            )
    else:
        files = [pattern]
    return files


def get_fastp_pipe_input(wildcards):
    return get_raw_reads(wildcards.sample, wildcards.unit, wildcards.fq)


def get_fastqc_input(wildcards):
    return get_raw_reads(wildcards.sample, wildcards.unit, wildcards.fq)[0]


def get_fastp_adapters(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]
    try:
        adapters = unit["adapters"]
        if isinstance(adapters, str):
            # Autotrimming is enabled by default.
            # Therefore no adapter parameter needs to be passed.
            if adapters == "auto_trim":
                return ""
            else:
                return adapters
        return ""
    except KeyError:
        return ""


def get_fastp_extra(wildcards):
    extra = config["params"]["fastp"]
    if "umi_read" in samples.columns and "umi_len" in samples.columns:
        if sample_has_umis(wildcards.sample):
            umi_extra = get_annotate_umis_params(wildcards)
            if extra and umi_extra:
                extra += " " + umi_extra
            else:
                extra += umi_extra
    return extra


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired


def get_map_reads_input(wildcards):
    if is_paired_end(wildcards.sample):
        return [
            f"results/merged/{wildcards.sample}_R1.fastq.gz",
            f"results/merged/{wildcards.sample}_R2.fastq.gz",
        ]
    return f"results/merged/{wildcards.sample}_single.fastq.gz"


def get_trimming_input(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    return f"results/mapped/vg/{wildcards.sample}.sorted.{ext}"


def get_primer_bed(wc):
    if isinstance(primer_panels, pd.DataFrame):
        if not pd.isna(primer_panels.loc[wc.panel, "fa2"]):
            return "results/primers/{}_primers.bedpe".format(wc.panel)
        else:
            return "results/primers/{}_primers.bed".format(wc.panel)
    else:
        if config["primers"]["trimming"].get("primers_fa2", ""):
            return "results/primers/uniform_primers.bedpe"
        else:
            return "results/primers/uniform_primers.bed"


def extract_unique_sample_column_value(sample, col_name):
    result = samples.loc[samples["sample_name"] == sample, col_name].drop_duplicates()
    if len(result) > 1:
        raise ValueError(
            "If a sample is specified multiple times in a samples.tsv "
            "sheet, all columns except 'group' must contain identical "
            "entries across the occurrences (rows).\n"
            f"Here we have sample '{sample}' with multiple entries for "
            f"the '{col_name}' column, namely:\n"
            f"{result}\n"
        )

    result = result.squeeze()
    return result


def get_sample_primer_fastas(sample):
    if isinstance(primer_panels, pd.DataFrame):
        panel = extract_unique_sample_column_value(sample, "panel")
        if not pd.isna(primer_panels.loc[panel, "fa2"]):
            return [
                primer_panels.loc[panel, "fa1"],
                primer_panels.loc[panel, "fa2"],
            ]
        return primer_panels.loc[panel, "fa1"]
    else:
        if config["primers"]["trimming"].get("primers_fa2", ""):
            return [
                config["primers"]["trimming"]["primers_fa1"],
                config["primers"]["trimming"]["primers_fa2"],
            ]
        return config["primers"]["trimming"]["primers_fa1"]


def get_panel_primer_input(panel):
    if panel == "uniform":
        if config["primers"]["trimming"].get("primers_fa2", ""):
            return [
                config["primers"]["trimming"]["primers_fa1"],
                config["primers"]["trimming"]["primers_fa2"],
            ]
        return config["primers"]["trimming"]["primers_fa1"]
    else:
        panel = primer_panels.loc[panel]
        if not pd.isna(panel["fa2"]):
            return [panel["fa1"], panel["fa2"]]
        return panel["fa1"]


def get_primer_regions(wc):
    if isinstance(primer_panels, pd.DataFrame):
        panel = extract_unique_sample_column_value(wc.sample, "panel")
        return f"results/primers/{panel}_primer_regions.tsv"
    return "results/primers/uniform_primer_regions.tsv"


def get_read_group(prefix: str):
    def inner(wildcards):
        """Denote sample name and platform in read group."""
        platform = extract_unique_sample_column_value(wildcards.sample, "platform")
        return "{prefix}'@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
            sample=wildcards.sample, platform=platform, prefix=prefix
        )

    return inner


def get_trimmed_fastqs(wc):
    adapters = units.loc[wc.sample, "adapters"]
    if adapters.notna().any() and adapters.isna().any():
        error_msg = (
            "Mixed adapter configuration for sample {!r}: units={!r}. "
            "Some units have adapters set and others are NA. "
            "Fix units.tsv so adapters is set for all units of the sample or NA for all."
        ).format(wc.sample, list(units.loc[wc.sample, "unit_name"]))
        raise Exception(error_msg)

    if adapters.notna().all():
        return expand(
            "results/trimmed/{sample}/{unit}_{read}.fastq.gz",
            unit=units.loc[wc.sample, "unit_name"],
            sample=wc.sample,
            read=wc.read,
        )
    else:
        fq = "fq1" if wc.read == "R1" or wc.read == "single" else "fq2"
        return [
            read
            for unit in units.loc[wc.sample, "unit_name"]
            for read in get_raw_reads(wc.sample, unit, fq)
        ]


def sample_has_umis(sample):
    return pd.notna(extract_unique_sample_column_value(sample, "umi_read"))


def get_annotate_umis_params(wildcards):
    translate_param = {"fq1": "read1", "fq2": "read2", "both": "per_read"}
    return "--umi --umi_loc {read} --umi_len {umi_len}".format(
        read=translate_param[
            extract_unique_sample_column_value(wildcards.sample, "umi_read")
        ],
        umi_len=str(extract_unique_sample_column_value(wildcards.sample, "umi_len")),
    )


def get_filter_params(wc):
    if isinstance(get_panel_primer_input(wc.panel), list):
        return "-b -F 12"
    return "-b -F 4"


def get_single_primer_flag(wc):
    if not isinstance(get_sample_primer_fastas(wc.sample), list):
        return "--first-of-pair"
    return ""


def get_shortest_primer_length(primers):
    primers = primers if isinstance(primers, list) else [primers]
    # set to 32 to match bwa-mem default value considering offset of 2
    min_length = 32
    for primer_file in primers:
        with open(primer_file, "r") as primer_f:
            lines = primer_f.readlines()
        min_primer = min(
            len(line.strip()) for i, line in enumerate(lines) if i % 2 == 1
        )
        min_length = min(min_length, min_primer)
    return min_length


def get_primer_extra(wc, input):
    extra = f"-R '@RG\tID:{wc.panel}\tSM:{wc.panel}' -L 100"
    min_primer_len = get_shortest_primer_length(input.reads)
    # Check if shortest primer is below default values
    if min_primer_len < 32:
        extra += f" -T {max(1, min_primer_len-2)}"
    if min_primer_len < 19:
        extra += f" -k {min_primer_len}"
    return extra


def get_fastqc_results(wildcards):
    group_samples = get_group_samples(wildcards.group)
    sample_units = units.loc[group_samples]
    sra_units = pd.isna(sample_units["fq1"])
    paired_end_units = sra_units | ~pd.isna(sample_units["fq2"])

    # fastqc
    pattern = "results/qc/fastqc/{unit.sample_name}/{unit.unit_name}.{fq}_fastqc.zip"
    yield from expand(pattern, unit=sample_units.itertuples(), fq="fq1")
    yield from expand(
        pattern, unit=sample_units[paired_end_units].itertuples(), fq="fq2"
    )

    # fastp
    if sample_units["adapters"].notna().all():
        pattern = "results/trimmed/{unit.sample_name}/{unit.unit_name}.{mode}.qc.html"
        yield from expand(
            pattern, unit=sample_units[paired_end_units].itertuples(), mode="paired"
        )
        yield from expand(
            pattern, unit=sample_units[~paired_end_units].itertuples(), mode="single"
        )

    # samtools idxstats (CRAM)
    yield from expand(
        "results/qc/{sample}.cram.idxstats",
        sample=group_samples,
    )

    # samtools stats (CRAM)
    yield from expand(
        "results/qc/{sample}.cram.stats",
        sample=group_samples,
    )


def get_pangenome_url(datatype):
    build = config["ref"]["build"].lower()
    source = config["ref"]["pangenome"]["source"]
    version = config["ref"]["pangenome"]["version"]
    if config["ref"]["species"] != "homo_sapiens" or build not in ["grch37", "grch38"]:
        raise ValueError(
            "Unsupported combination of species and build. Only homo_sapiens and GRCh37/GRCh38 are supported for pangenome mapping."
        )
    if source != "hprc":
        raise ValueError(
            "Unsupported pangenome source. Only 'hprc' is currently supported."
        )

    base_url = config["ref"]["pangenome"].get(
        "base_url",
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus",
    )

    return f"{base_url}/hprc-{version}-mc-{build}/hprc-{version}-mc-{build}.{datatype}"
