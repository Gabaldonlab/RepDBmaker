# RepDBmaker Pipeline


<p align="center">
  <img src="images/Fig1.png" alt="pipeline_schema" width="800"/>
</p>

This repository contains a Snakemake-based workflow to build taxonomically annotated protein sequence databases from multiple public sources and custom genome collections.

## Table of contents

- [RepDBmaker Pipeline](#repdbmaker-pipeline)
  - [Table of contents](#table-of-contents)
  - [What this pipeline does](#what-this-pipeline-does)
  - [Repository layout](#repository-layout)
  - [Requirements](#requirements)
  - [Quick start](#quick-start)
  - [Configuration](#configuration)
  - [Pipeline components](#pipeline-components)
    - [Downloads and metadata](#downloads-and-metadata)
    - [Taxonomy assembly](#taxonomy-assembly)
    - [Clustering](#clustering)
    - [Decontamination](#decontamination)
  - [Outputs](#outputs)
  - [Benchmarking](#benchmarking)
  - [Custom databases](#custom-databases)
  - [Useful files](#useful-files)
  - [Authors](#authors)
  - [References](#references)
  - [TODO](#todo)

## What this pipeline does

- Downloads and assembles taxonomic resources for:
  - prokaryotes from [GTDB](https://gtdb.ecogenomic.org/)
  - eukaryotes from [EukProt](https://evocellbio.com/eukprot/), [P10K](https://ngdc.cncb.ac.cn/p10k/) and [UniProt](https://www.uniprot.org/proteomes?query=%28taxonomy_id%3A2759%29)
  - viruses from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/)
- Parses genomes into protein FASTA files with accession-to-taxonomy mapping
- Creates searchable sequence databases in Diamond, MMseqs2, and BLAST formats
- Optionally clusters database sequences by taxonomic rank
- Optionally decontaminates selected databases using sequence cluster analysis

## Repository layout

- `config/repdb.yaml` - main pipeline configuration
- `workflow/Snakefile` - main Snakemake entry point
- `workflow/rules/` - rule files for download, taxonomy, database building, clustering, decontamination, and stats
- `workflow/envs/` - Conda environment definitions used by the workflow
- `resources/` - static input files and genome lists used by the pipeline

## Requirements

- Snakemake
- Conda or Miniconda
- Internet access for downloading external resources
- Disk space for large genomic and database collections

The workflow uses the following Conda environments:

- `workflow/envs/python.yaml` — Python, pandas, matplotlib, polars
- `workflow/envs/homology.yaml` — diamond, mmseqs2, blast
- `workflow/envs/utils.yaml` — taxonkit, csvtk, ncbi-datasets-cli, jq, seqkit
- `workflow/envs/R.yaml` — R and visualization/taxonomy packages

## Quick start

1. Load Conda/Snakemake and activate your environment:

```bash
conda activate snakemake
```

2. In case you are interested in the available proteomes you can just un this:

```bash
snakemake -j 1 --until available_proteomes
```

It will produce the `results/meta/available_proteomes.tsv` that can be easily parsed to further select whatever susbet of proteomes you may be interested in.

3. Alternatively you can run the full workflow directly:

```bash
snakemake -j 14
```

Sometimes things can go wrong while downloading a proteome. That is why there is a step (rule `db_stats`) that will fail if any gzipped fasta is malformed and will block the creation of the database fasta. 

In this case you can run for example:

```bash
cut -f1 results/dbs/repdb/genome_table.tsv | xargs -I {} sh -c 'gzip -t "{}" || echo "Failed: {}"'
```

Then delete the problematic ones and re-run the pipeline. If the problem persitsts, there may be other sort of problems
(this specific proteome files are broken or the current downloading script fails), I reccomend to exclude them and find the most suitable alternative. 


## Configuration

Edit `config/repdb.yaml` to select which databases and formats to build. You can also use it as template for any `yaml` file as long as you run:

`snakemake --configfile path/to/custom.yaml`

Key configuration sections:

- `dbs.type` — allowed values: `diamond`, `mmseqs`, `blastp`
- `dbs.build.repdb` — enable RepDB construction
- `dbs.build.custom` — define custom databases with genome IDs and optional clustering/decontamination
- `dbs.gtdb_version` — GTDB release version, `default` is the latest.
- `files.clades_to_keep` — eukaryotic lineages to retain
- `files.genomes_to_exclude` — genomes excluded from RepDB
- `files.new_genomes` — custom genome metadata for `results/proteomes/cus/`

Example custom database entry:

```yaml
dbs:
  type: ["diamond", "mmseqs", "blastp"]
  build:
    smalleuks:
      ids: resources/eukas_50.ids
      cluster:
        level: order    # Taxonomic Rank to constrain clustering 
        identity: 0.8   # Clustering Identity
        coverage: 0.8   # Clustering coverage
      decontaminate:
        identity: 0.9   # Clustering identity
        coverage: 0.5   # Clustering coverage
        cov_mode: 3     # Clustering coverage mode
        prop_euka: 0.5  # % of Eukas in contaminant clusters
```

## Pipeline components

### Downloads and metadata

The workflow retrieves:

- NCBI taxdump
- GTDB taxonomy and representative genomes
- NCBI RefSeq viral genomes
- UniProt reference proteome metadata
- EukProt metadata and proteomes
- P10K metadata and lineage data
- UniEuk taxonomy resource

### Taxonomy assembly

- Builds combined taxonomic annotations for viruses, GTDB, and selected eukaryotes
- Creates RepDB-specific taxonomy using filtered eukaryotes, all viruses, and all GTDB genomes
- Generates a full taxdump used by database builders and clustering steps


### Clustering

- Optionally clusters database sequences by taxonomic clade
- Uses taxonomic rank values such as `class`, `order`, `family`, etc.
- Generates per-clade FASTA and cluster output, then merges into a clustered database file

### Decontamination

- Concatenates RepDB and custom database FASTA files for decontamination when configured
- Runs MMseqs2 clustering to identify sequence clusters
- Extracts non-singleton clusters and mixed clusters with eukaryotic content
- Produces `contaminants.txt`, `pair_counts.tsv`, and other diagnostic outputs

## Outputs

Primary output directories and files:

- `results/dbs/<db>/` — built database outputs
  - `<db>.fa.gz` — compressed protein FASTA
  - `<db>_accession_map.txt` — accession-to-taxid mapping
  - `<db>_map` and `<db>_nohead.map` — BLAST/MMseqs maps
  - `<db>_diamond` — Diamond database index
  - `<db>_mmseqs` — MMseqs2 database index
  - `<db>_blastp` — BLASTP database files
- `results/dbs/<db>/cluster/cluster_params.yaml` — clustering settings
- `results/dbs/<db>/cluster/<db>_clustered.fa.gz` — clustered FASTA output
- `results/dbs/<db>/decontaminate/decontaminate_params.yaml` — decontamination settings
- `results/dbs/<db>/decontaminate/contaminants.txt` — decontamination candidates
- `results/dbs/<db>/decontaminate/pair_counts.tsv` — cluster pair counts
- `results/stats/` — database and clustering statistics
- `results/meta/check_resources.txt` — validation of downloaded assets
- `results/taxonomies/` — taxonomy annotations and RepDB selection outputs
- `results/taxdump/repdb_taxdump/` — taxonkit taxdump for combined genomes

## Benchmarking

A separate benchmark workflow is available at `workflow/benchmark.smk` and uses `config/benchmark.yaml` to compare RepDB against reference databases such as NR and clustered NR.

## Custom databases

To add a custom database:

1. Define an entry under `dbs.build.custom` in `config/repdb.yaml`
2. Provide an `ids` file containing the genome identifiers or metadata
3. Optionally add `cluster` and `decontaminate` sections for clustering and contamination filtering

The pipeline will create `results/dbs/<custom_db>/` with the expected outputs.


## Useful files

- `workflow/Snakefile` — main workflow logic
- `workflow/rules/` — modular rule definitions
- `workflow/envs/` — Conda environments used by rules
- `config/repdb.yaml` — pipeline configuration
- `resources/` — genome lists, excluded genomes, and clade filters


## Authors

- Firstname Lastname
  - Affiliation
  - ORCID profile
  - home page

## References

> ADD

## TODO

- UniProt: include also non kb proteomes?
- JGI?
- unit tests
- proper README
- docker container
- snakevision
- images folder
- profiles
- schemas
- github workflows (check snakemake template)
