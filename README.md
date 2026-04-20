# RepDBmaker Pipeline

<p align="center">
  <img src="images/Fig1.png" alt="pipeline_schema" width="800"/>
</p>

A Snakemake workflow for building taxonomically annotated protein sequence databases from public resources and custom genome collections.

## Overview

RepDBmaker assembles protein sequence databases from multiple sources, annotates them with taxonomy, and builds searchable indices for:

- Diamond
- MMseqs2
- BLAST

It supports:

- prokaryotes from [GTDB](https://gtdb.ecogenomic.org/)
- eukaryotes from [EukProt](https://evocellbio.com/eukprot/), [P10K](https://ngdc.cncb.ac.cn/p10k/), and [UniProt](https://www.uniprot.org/proteomes?query=%28taxonomy_id%3A2759%29)
- viruses from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/)

Optional features include taxonomic clustering and contamination filtering.


## Insallation

### Requirements

- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- Conda or Miniconda
- Internet access for external downloads
- Sufficient disk space for genome and database files

If using `--sdm conda`, Snakemake will automatically create the following environments from `workflow/envs/`:

- `workflow/envs/python.yaml` — Python, pandas, matplotlib, polars
- `workflow/envs/homology.yaml` — diamond, mmseqs2, blast
- `workflow/envs/utils.yaml` — taxonkit, csvtk, ncbi-datasets-cli, jq, seqkit
- `workflow/envs/R.yaml` — R and visualization/taxonomy packages

Create all environments without executing the pipeline:

```bash
snakemake --conda-create-envs-only
```

### Docker

A Docker image is available at [Docker Hub](https://hub.docker.com/repository/docker/gmuttiirb/repdbmaker/general).

If [Docker](https://docs.docker.com/engine/install/) is installed, run the pipeline from the repository root with:

```bash
docker run --rm -v $(pwd):/app/data gmuttiirb/repdbmaker:v1.0 snakemake --cores 2 --directory /app/data -n
```

## Quick start

To inspect the available proteomes before running the full workflow:

```bash
snakemake -j 1 --until available_proteomes
```

This generates `results/meta/available_proteomes.tsv`, which can be used to select a proteome subset.

To run the full workflow:

```bash
snakemake -j 14
```

## Troubleshooting

Sometimes things can go wrong while downloading a proteome. That is why there is a step (rule `db_stats`) that will fail if any gzipped fasta is malformed and will block the creation of the database fasta. 

To check for broken files:

```bash
cut -f1 results/dbs/<db>/genome_table.tsv | xargs -I {} sh -c 'gzip -t "{}" || echo "Failed: {}"'
```

Then delete the problematic ones and re-run the pipeline. If the problem persists, there may be other sort of problems
(the files may be broken or the current downloading script fails), I reccomend to exclude them and find the most suitable alternative. 

## Configuration

Configure the workflow in `config/repdb.yaml` or pass a custom config file with:

```bash
snakemake --configfile path/to/custom.yaml
```

Example configuration:

```yaml
# Database configuration
dbs:
  gtdb_version: "latest"  # use release220, release226, etc.
  type: ["diamond", "mmseqs"]  # Available types: diamond, mmseqs, blastp
  build:
    repdb:
      decontamination:
        # optimized settings for ContScout benchmarking
        identity: 0.9
        coverage: 0.5
        cov_mode: 3
        prop_euka: 0.5
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

Using this config file will allow the creation of RepDB and its clustered version.

## Custom databases

To add any custom database:

1. Define an entry under `dbs.build.custom`
2. Provide an `ids` file with genome identifiers or metadata
3. Optionally configure `cluster` and `decontaminate`

Example:

```yaml
dbs:
  type: ["diamond", "mmseqs", "blastp"]
  build:
    smalleuks:
      ids: resources/eukas_50.ids
      cluster:
        level: order
        identity: 0.8
        coverage: 0.8
      decontaminate:
        identity: 0.9
        coverage: 0.5
        cov_mode: 3
        prop_euka: 0.5
```

The workflow will create `results/dbs/<custom_db>/` and its associated outputs.

## Outputs

Key output locations:

- `results/dbs/<db>/`
  - `<db>.fa.gz` — compressed protein FASTA
  - `<db>_accession_map.txt` — accession-to-taxid mapping
  - `<db>_map`, `<db>_nohead.map` — BLAST/MMseqs maps
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
- `results/taxonomies/` — taxonomy annotations and selection outputs
- `results/taxdump/repdb_taxdump/` — taxonkit taxdump for combined genomes

## Utilities

The repository includes helper scripts for working with generated databases.

- `workflow/scripts/get_fasta.py`: extract a subset of FASTA sequences from an MMseqs2 database

Example usage:

```bash
python workflow/scripts/get_fasta.py <id_file> <output_fasta> <db_mmseqs>
```

## Benchmarking

A benchmark workflow is available in `workflow/benchmark.smk` and uses `config/benchmark.yaml` to compare RepDB against reference databases such as NR and clustered NR.

A results notebook is available at `workflow/notebooks/comparison.Rmd`.

## Citation

Please add the preferred citation here.

## License

See the `LICENSE` file for license details.
