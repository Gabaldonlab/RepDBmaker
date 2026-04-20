# RepDBmaker Pipeline


<p align="center">
  <img src="images/Fig1.png" alt="pipeline_schema" width="800"/>
</p>

This repository contains a Snakemake-based workflow to build taxonomically annotated protein sequence databases from multiple public sources and custom genome collections.

## The pipeline

- Downloads and assembles taxonomic resources for:
  - prokaryotes from [GTDB](https://gtdb.ecogenomic.org/)
  - eukaryotes from [EukProt](https://evocellbio.com/eukprot/), [P10K](https://ngdc.cncb.ac.cn/p10k/) and [UniProt](https://www.uniprot.org/proteomes?query=%28taxonomy_id%3A2759%29)
  - viruses from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/)
- Parses genomes into protein FASTA files with accession-to-taxonomy mapping
- Creates searchable sequence databases in Diamond, MMseqs2, and BLAST formats
- Optionally clusters database sequences by taxonomic rank
- Optionally decontaminates selected databases using sequence cluster analysis

## Installation

### Reuirements

- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- Conda or Miniconda
- Internet access for downloading external resources
- Disk space for large genomic and database collections

The workflow will then automatically download the following Conda environments (if `--sdm conda` is used):

- `workflow/envs/python.yaml` — Python, pandas, matplotlib, polars
- `workflow/envs/homology.yaml` — diamond, mmseqs2, blast
- `workflow/envs/utils.yaml` — taxonkit, csvtk, ncbi-datasets-cli, jq, seqkit
- `workflow/envs/R.yaml` — R and visualization/taxonomy packages

You can run this command to create all the necessary environments without running the pipeline:

```bash
snakemake --conda-create-envs-only 
```

### Docker image

Alternatively, a Docker image is available at [Docker Hub](https://hub.docker.com/repository/docker/gmuttiirb/repdbmaker/general).

In this case the only dependency will be installing [Docker](https://docs.docker.com/engine/install/).

The pipeline can be used with this command: 

```bash
docker run --rm -v $(pwd):/app/data gmuttiirb/repdbmaker:v1.0 snakemake --cores 2 --directory /app/data -n
```

## Getting started


If first, you are interested in the available proteomes you can just run this:

```bash
snakemake -j 1 --until available_proteomes
```

It will produce the `results/meta/available_proteomes.tsv` that can be easily parsed to further select whatever susbet of proteomes you may be interested in.

Alternatively you can run the full workflow directly:

```bash
snakemake -j 14
```

Sometimes things can go wrong while downloading a proteome. That is why there is a step (rule `db_stats`) that will fail if any gzipped fasta is malformed and will block the creation of the database fasta. 

If there are any issues, you can run for example:

```bash
cut -f1 results/dbs/repdb/genome_table.tsv | xargs -I {} sh -c 'gzip -t "{}" || echo "Failed: {}"'
```

Then delete the problematic ones and re-run the pipeline. If the problem persists, there may be other sort of problems
(the files may be broken or the current downloading script fails), I reccomend to exclude them and find the most suitable alternative. 

## Configuration

Edit `config/repdb.yaml` to select which databases and formats to build. You can also use it as template for any `yaml` file as long as you run:

`snakemake --configfile path/to/custom.yaml`

The default configuration file looks like this:

```yaml
# Database configuration
dbs:
  gtdb_version: "latest"  # use release220 or release226 etc to use older versions
  type: ["diamond", "mmseqs"]  # Available database types: diamond, mmseqs, blastp
  # Which databases to build
  build:
    repdb: 
      decontamination:
        # these parameters were optimized through a benchmark with ContScout 
        identity: 0.9
        coverage: 0.5
        cov_mode: 3
        prop_euka: 0.5
    # Custom databases: uncomment and add your database entries
    # Example:
    custom:
      clusteredrepdb: 
        ids: resources/repdb.ids
        cluster: 
          level: class
          identity: 0.9
          coverage: 0.9

files: 
  clades_to_keep: "resources/clades_tokeep.txt"
  genomes_to_exclude: "resources/exclude.txt"
  new_genomes: "resources/custom_genomes_repdb.csv"
```

It will allow the creation of RepDB and its clustered version.

### Custom databases

To add any custom database:

1. Define an entry under `dbs.build.custom` in the `config` file
2. Provide an `ids` file containing the genome identifiers or metadata
3. Optionally add `cluster` and `decontaminate` sections for clustering and contamination filtering

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

The pipeline will create `results/dbs/<custom_db>/` with the expected outputs.

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

A separate benchmark workflow is available at `workflow/benchmark.smk` and uses `config/benchmark.yaml` to compare RepDB against reference databases such as NR and clustered NR. A R Markdown notebook showing our results is available in `workflow/notebooks/comparison.Rmd`.

## Authors

- Firstname Lastname
  - Affiliation
  - ORCID profile
  - home page

## Citation

> ADD

## TODO

- UniProt: include also non kb proteomes?
- JGI?
- unit tests
- proper README
- images folder
- profiles
- schemas
- github workflows (check snakemake template)
