import glob

import pandas as pd
from snakemake.utils import validate
from snakemake.exceptions import WorkflowError

validate(config, schema="../schemas/config.schema.yaml")


samples = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype={"sample_name": str, "umi_len": "Int64"},
        comment="#",
    )
    .set_index("sample_name", drop=False)
    .sort_index()
)


wildcard_constraints:
    sample="|".join(samples["sample_name"]),


# construct genome name
datatype_genome = "dna"
species = config["ref"]["species"]
build = config["ref"]["build"]
release = config["ref"]["release"]
genome_name = f"genome.{datatype_genome}.{species}.{build}.{release}"
genome_prefix = f"<resources>/{genome_name}"
genome = f"{genome_prefix}.fasta"
genome_fai = f"{genome}.fai"
genome_dict = f"{genome_prefix}.dict"
pangenome_name = f"pangenome.{species}.{build}"
pangenome_prefix = f"<resources>/{pangenome_name}"


def _group_or_sample(row):
    group = row.get("group", None)
    if pd.isnull(group):
        return row["sample_name"]
    return group


samples["group"] = [_group_or_sample(row) for _, row in samples.iterrows()]

if "umi_read" not in samples.columns:
    samples["umi_read"] = pd.NA

validate(samples, schema="../schemas/samples.schema.yaml")


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
    (
        pd.read_csv(
            config["primers"]["trimming"]["tsv"],
            sep="\t",
            dtype={"panel": str, "fa1": str, "fa2": str},
            comment="#",
        )
        .set_index(["panel"], drop=False)
        .sort_index()
    )
    if config["primers"]["trimming"].get("tsv", "")
    else None
)


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def get_results_bqsr():
    if lookup(
        within=config,
        dpath="base_recalibration/activate",
        default=False,
    ):
        return "results"
    return "results/recal"


def get_results_consensus():
    if (
        lookup(
            within=config,
            dpath="calc_consensus_reads/activate",
            default=False,
        )
        and not get_results_bqsr()
    ):
        return "results"
    return "results/consensus"


def get_results_trimmed():
    primers_active = (
        bool(config["primers"]["trimming"].get("primers_fa1"))
        or bool(config["primers"]["trimming"].get("primers_fa2"))
        or ("panel" in samples.columns and samples["panel"].notna().any())
    )

    if (
        primers_active
        and not bool(
            lookup(
                within=config,
                dpath="calc_consensus_reads/activate",
                default=False,
            )
        )
        and not bool(
            lookup(
                within=config,
                dpath="base_recalibration/activate",
                default=False,
            )
        )
    ):
        return "results"
    return "results/trimmed"


def get_results_dedup():
    primers_active = (
        bool(config["primers"]["trimming"].get("primers_fa1"))
        or bool(config["primers"]["trimming"].get("primers_fa2"))
        or ("panel" in samples.columns and samples["panel"].notna().any())
    )

    if (
        bool(
            lookup(
                within=config,
                dpath="remove_duplicates/activate",
                default=False,
            )
        )
        and not primers_active
        and not bool(
            lookup(
                within=config,
                dpath="calc_consensus_reads/activate",
                default=False,
            )
        )
        and not bool(
            lookup(
                within=config,
                dpath="base_recalibration/activate",
                default=False,
            )
        )
    ):
        return "results"
    return "results/dedup"


def mapped_stage_is_final():
    primers_active = (
        bool(config["primers"]["trimming"].get("primers_fa1"))
        or bool(config["primers"]["trimming"].get("primers_fa2"))
        or ("panel" in samples.columns and samples["panel"].notna().any())
    )

    return (
        not bool(
            lookup(
                within=config,
                dpath="remove_duplicates/activate",
                default=False,
            )
        )
        and not primers_active
        and not bool(
            lookup(
                within=config,
                dpath="calc_consensus_reads/activate",
                default=False,
            )
        )
        and not bool(
            lookup(
                within=config,
                dpath="base_recalibration/activate",
                default=False,
            )
        )
    )


def get_results_mapped_bwa():
    if not is_activated("ref/pangenome") and mapped_stage_is_final():
        return "results"
    return "results/mapped/bwa"


def get_results_mapped_vg():
    if is_activated("ref/pangenome") and mapped_stage_is_final():
        return "results"
    return "results/mapped/vg"


get_aligner = branch(
    is_activated("ref/pangenome"),
    then="vg",
    otherwise="bwa",
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
            return f"pipe/fastp/{unit.sample_name}/{unit.unit_name}.{fq}.fastq{ending}"

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
    return expand("sra/{accession}_{read}.fastq.gz", accession=accession, read=fq)


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


get_map_reads_input = branch(
    lambda wc: is_paired_end(wc.sample),
    then=[
        "<results>/merged/{sample}_R1.fastq.gz",
        "<results>/merged/{sample}_R2.fastq.gz",
    ],
    otherwise="<results>/merged/{sample}_single.fastq.gz",
)


def get_mapped_stage_input(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    if is_activated("ref/pangenome"):
        return f"<results_mapped_vg>/{{sample}}.{ext}"
    return f"<results_mapped_bwa>/{{sample}}.{ext}"


def get_recalibrate_quality_input(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    if is_activated("calc_consensus_reads"):
        return f"<results_consensus>/{{sample}}.{ext}"
    return get_consensus_input(wildcards, bai)


def get_consensus_input(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    if sample_has_primers(wildcards):
        return f"<results_trimmed>/{{sample}}.{ext}"
    return get_trimming_input(wildcards, bai)


def get_trimming_input(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    if is_activated("remove_duplicates"):
        return f"<results_dedup>/{{sample}}.{ext}"
    return get_mapped_stage_input(wildcards, bai)


def extract_unique_sample_column_value(sample, col_name):
    result = samples.loc[samples["sample_name"] == sample, col_name].drop_duplicates()
    if type(result) is not str:
        if len(result) > 1:
            raise ValueError(
                "If a sample is specified multiple times in a samples.tsv"
                "sheet, all columns except 'group' must contain identical"
                "entries across the occurrences (rows).\n"
                f"Here we have sample '{sample}' with multiple entries for"
                f"the '{col_name}' column, namely:\n"
                f"{result}\n"
            )
        else:
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


get_primer_regions = branch(
    isinstance(primer_panels, pd.DataFrame),
    then=lambda wc: f"<results>/primers/{extract_unique_sample_column_value(wc.sample, 'panel')}_primer_regions.tsv",
    otherwise="<results>/primers/uniform_primer_regions.tsv",
)


def get_markduplicates_extra(wildcards):
    c = config["params"]["picard"]["MarkDuplicates"]

    if sample_has_umis(wildcards):
        b = "--BARCODE_TAG BX"
    else:
        b = ""

    if is_activated("calc_consensus_reads"):
        d = "--TAG_DUPLICATE_SET_MEMBERS true"
    else:
        d = ""

    return f"{c} {b} {d}"


def get_read_group(prefix: str):
    def inner(wildcards):
        """Denote sample name and platform in read group."""
        platform = extract_unique_sample_column_value(wildcards.sample, "platform")
        return r"{prefix}'@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
            sample=wildcards.sample, platform=platform, prefix=prefix
        )

    return inner


def sample_has_primers(wildcards):
    sample_name = wildcards.sample

    if config["primers"]["trimming"].get("primers_fa1") or (
        "panel" in samples.columns
        and samples.loc[samples["sample_name"] == sample_name, "panel"].notna().any()
    ):
        if not is_paired_end(sample_name):
            raise WorkflowError(
                f"Primer trimming is only available for paired-end data. Sample '{sample_name}' is not paired-end."
            )
        return True
    return False


def sample_has_umis(wildcards):
    return pd.notna(extract_unique_sample_column_value(wildcards.sample, "umi_read"))


get_fastp_extra = branch(
    sample_has_umis,
    then=lambda wc: config["params"]["fastp"] + get_annotate_umis_params(wc),
    otherwise=lambda wc: config["params"]["fastp"],
)


def get_annotate_umis_params(wildcards):
    translate_param = {"fq1": "read1", "fq2": "read2", "both": "per_read"}
    return "--umi --umi_loc {read} --umi_len {umi_len}".format(
        read=translate_param[
            extract_unique_sample_column_value(wildcards.sample, "umi_read")
        ],
        umi_len=str(extract_unique_sample_column_value(wildcards.sample, "umi_len")),
    )


def get_shortest_primer_length(primers):
    primers = primers if isinstance(primers, list) else [primers]
    # set to 32 to match bwa-mem default value considering offset of 2
    min_length = 32
    for primer_file in primers:
        with open(primer_file, "r") as p:
            min_primer = min(
                [len(p.strip()) for i, p in enumerate(p.readlines()) if i % 2 == 1]
            )
            min_length = min(min_length, min_primer)
    return min_length


def get_primer_extra(wildcards, input):
    extra = rf"-R '@RG\tID:{wildcards.panel}\tSM:{wildcards.panel}' -L 100"
    min_primer_len = get_shortest_primer_length(input.reads)
    # Check if shortest primer is below default values
    if min_primer_len < 32:
        extra += f" -T {min_primer_len-2}"
    if min_primer_len < 19:
        extra += f" -k {min_primer_len}"
    return extra


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
    return (
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/"
        "freeze1/minigraph-cactus/"
        f"hprc-{version}-mc-{build}/hprc-{version}-mc-{build}.{datatype}"
    )
