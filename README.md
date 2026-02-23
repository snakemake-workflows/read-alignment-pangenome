# Snakemake workflow: read-alignment-pangenome

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/snakemake-workflows/read-alignment-pangenome/workflows/Tests/badge.svg?branch=main)](https://github.com/snakemake-workflows/read-alignment-pangenome/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/snakemake-workflows/read-alignment-pangenome)

A Snakemake workflow for aligning sequencing reads to a pangenome graph with **vg giraffe**, producing **sorted CRAM + CRAI** per sample.


- [Snakemake workflow: read-alignment-pangenome](#snakemake-workflow-read-alignment-pangenome)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [Workflow profiles](#workflow-profiles)
  - [Authors](#authors)
  - [References](#references)
  - [TODO](#todo)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/snakemake-workflows/read-alignment-pangenome).

Detailed information about input data and workflow configuration can also be found in [`config/README.md`](config/README.md).

**Main target and outputs**

The main alignment target is:

- `only_alignment`

Expected outputs (per sample):

- `<results>/mapped/vg/{sample}.sorted.cram`
- `<results>/mapped/vg/{sample}.sorted.cram.crai`

> Note: If `adapters` is empty/NA in `units.tsv`, trimming with fastp is bypassed and raw reads are used.

## Deployment options

To run the workflow from command line, change the working directory.

```bash
cd path/to/read-alignment-pangenome
```

Adjust options in the default config file `config/config.yaml`.
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake -n --cores 1 --use-conda only_alignment
```

To run the workflow with **conda**:

```bash
snakemake --cores 2 --use-conda only_alignment
```

To run the workflow with **apptainer** / **singularity** (optional), if containers are configured (e.g., via a workflow profile and/or rule-level `container:` directives):

```bash
snakemake --cores 2 --use-conda --use-apptainer only_alignment
```

## Workflow profiles

The `profiles/` directory can contain any number of [workflow-specific profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) that users can choose from.
The [profiles `README.md`](profiles/README.md) provides more details.

## Authors

- Firstname Lastname
  - Affiliation
  - ORCID profile
  - home page

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

## TODO


- Update the [deployment](#deployment-options), [authors](#authors) and [references](#references) sections.
- Update the `README.md` badges. Add or remove badges for `conda`/`singularity`/`apptainer` usage depending on the workflow's [deployment](#deployment-options) options.

