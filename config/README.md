# General settings

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.

This workflow is derived from `snakemake-workflows/dna-seq-varlociraptor` and focuses on:
- reference + optional pangenome resource preparation
- read preprocessing/merging
- read mapping, including optional pangenome-graph alignment with `vg giraffe`
- alignment postprocessing to produce final BAM + BAI

Variant calling, annotation, filtering, and reporting are intentionally out of scope.

# Sample sheet

Add samples to `config/samples.tsv`. For each sample, the columns `sample_name`, `platform`, `group`, and `datatype` must be defined.

- Samples within the same `group` can be treated as belonging together for aggregation logic that is retained from upstream.
- The `platform` column needs to contain the used sequencing platform (one of `CAPILLARY`, `LS454`, `ILLUMINA`, `SOLID`, `HELICOS`, `IONTORRENT`, `ONT`, `PACBIO`). This is required because the workflow adds read groups during alignment postprocessing.
- The `datatype` column is used by upstream-derived helper logic to determine the alignment branch and related processing.
- Optionally, a `panel` column can be provided. This is only relevant if primer trimming is enabled panel-wise (see [Primer trimming](#primer-trimming)); the value links a sample to a primer panel definition.
- Optionally, the columns `umi_read` and `umi_len` can be provided to enable UMI annotation (see [Annotating UMIs](#annotating-umis)).
  - `umi_read` can be `fq1`, `fq2`, or `both`.
  - `umi_len` is the number of bases (UMI length) to be annotated as UMI.

Missing values can be specified by empty columns or by writing `NA`. Lines can be commented out with `#`.

# Unit sheet

For each sample, add one or more sequencing units (runs, lanes, or replicates) to the unit sheet `config/units.tsv`.

- Each unit has a `unit_name` (lane/run/replicate ID).
- Each unit has a `sample_name`, which associates it with the biological sample it comes from.
  This information is used to merge all units of a sample before read mapping.

For each unit, you need to specify one of these input modes:

- `fq1` only for single-end reads (path to a FASTQ file)
- `fq1` and `fq2` for paired-end reads (paths to FASTQ files)
- `sra` only: specify an SRA accession (for example, `SRR...`). The workflow will download paired-end reads from SRA.

If both local files (`fq1`, `fq2`) and an SRA accession (`sra`) are available, the local files will be used.

## Adapters / trimming behavior (fastp)

Adapters can be configured in the `adapters` column by putting [fastp arguments](https://github.com/OpenGene/fastp?tab=readme-ov-file#adapters) in quotation marks
(for example, `"--adapter_sequence ACGC... --adapter_sequence_r2 GCTA..."`).

Automatic adapter trimming can be enabled by setting the keyword `auto_trim`.
If the `adapters` column is empty or `NA` for any unit of a sample, fastp will not be used for that sample and raw reads will be merged directly.

Missing values can be specified by empty columns or by writing `NA`. Lines can be commented out with `#`.

# Primer trimming

Primer trimming is retained from upstream logic. Trimming will be applied if global primer sequences are provided in `config/config.yaml` or primer panels are set in the sample sheet (column `panel`).

Primers can be defined either:

- directly in `config/config.yaml` (`primers.trimming.primers_fa1` / `primers.trimming.primers_fa2`), or
- via a separate TSV file (`primers.trimming.tsv`) with columns: `panel`, `fa1`, `fa2` (optional).

If a panel is not provided for a sample, primer trimming will not be performed for that sample.

For single-primer trimming only, only the first entry (`fa1`) needs to be defined.

# Annotating UMIs

UMI annotation is retained from upstream logic.

To enable it, add the following columns to `config/samples.tsv`:

- `umi_read`: where the UMI is located
  - `fq1` if the UMI is part of read 1
  - `fq2` if the UMI is part of read 2
  - `both` if there are UMIs in both paired-end reads
- `umi_len`: number of bases (UMI length) to be annotated as UMI.