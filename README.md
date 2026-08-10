# RepDBmaker Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A58.11.6-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Docker](https://img.shields.io/badge/docker-gmuttiirb%2Frepdbmaker-blue?logo=docker)](https://hub.docker.com/repository/docker/gmuttiirb/repdbmaker/general)

A Snakemake workflow for building taxonomically annotated protein sequence
databases from public resources and custom genome collections.

<p align="center">
  <img src="images/Fig1.png" alt="pipeline_schema" width="800"/>
</p>

## Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Reproducibility](#reproducibility)
- [Troubleshooting](#troubleshooting)
- [Documentation](#documentation)
- [Citation](#citation)

## Overview

RepDBmaker assembles protein sequence databases from multiple sources, annotates
them with taxonomy, and builds searchable indices for Diamond, MMseqs2 and BLAST.

It draws on:

- prokaryotes from [GTDB](https://gtdb.ecogenomic.org/)
- eukaryotes from [EukProt](https://evocellbio.com/eukprot/), [P10K](https://ngdc.cncb.ac.cn/p10k/), and [UniProt](https://www.uniprot.org/proteomes?query=%28taxonomy_id%3A2759%29)
- viruses from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/)
- your own genomes or transcriptomes, alongside the public sources

with optional taxonomic clustering and contamination filtering along the way.

## Installation

### Requirements

- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (v8.11.6)
- Conda or Miniconda
- Internet access for external downloads
- Sufficient disk space for genome and database files

Pick how the workflow schedules its jobs — local machine, SLURM, LSF, or a
site-specific profile — in [Choosing an executor](docs/executors.md).

If using `--sdm conda`, Snakemake will automatically create the following
environments from `workflow/envs/`:

- `workflow/envs/python.yaml` — Python, pandas, matplotlib, polars
- `workflow/envs/homology.yaml` — diamond, mmseqs2, blast
- `workflow/envs/utils.yaml` — taxonkit, csvtk, ncbi-datasets-cli, jq, seqkit
- `workflow/envs/R.yaml` — R and visualization/taxonomy packages

Create all of them up front, without running the pipeline:

```bash
snakemake --sdm conda --conda-create-envs-only
```

These `.yaml` specs are intentionally loose (minimum bounds only where a feature
requires it) so a fresh install resolves against current packages. To instead
reproduce the **exact** package builds used for the published RepDB v1.0, see
[Reproducibility](#reproducibility) — a pin file sitting next to each `.yaml` is
picked up automatically by `--sdm conda`, no extra flag needed.

### Docker

A Docker image is available at [Docker Hub](https://hub.docker.com/repository/docker/gmuttiirb/repdbmaker/general).
It's built from the same pin files, so it reproduces the RepDB v1.0 toolchain
rather than re-solving at build time. From the repository root:

```bash
docker run --rm -t -v $(pwd):/app/data gmuttiirb/repdbmaker:v1.0 \
  snakemake --configfile config/repdb.yaml --cores 2 --directory /app/data -n
```

For reproducible pulls, reference the image by its immutable digest instead of
the mutable `:v1.0` tag (`docker inspect --format='{{index .RepoDigests 0}}'
gmuttiirb/repdbmaker:v1.0` prints the current one):

```bash
docker run --rm -t -v $(pwd):/app/data \
  gmuttiirb/repdbmaker@sha256:1c0d3f397ac79ff48f712448156d1a60eb751290fb7becb446c12a11d3a791bc \
  snakemake --configfile config/repdb.yaml --cores 2 --directory /app/data -n
```

For an actual (non dry) run, add `--sdm conda --conda-prefix /conda-envs` —
both matter. Without `--sdm conda`, each rule's `conda:` directive is ignored
and the command runs against the base image's bare `PATH`, failing on the
first missing tool. `--conda-prefix /conda-envs` must be exactly this value:
it's not just a path setting, it also feeds the hash Snakemake uses to
recognize the image's pre-built envs — get it wrong and Snakemake rebuilds
everything from scratch under `.snakemake/conda/` in the bind-mounted,
empty `/app/data` instead.

```bash
docker run --rm -t -v $(pwd):/app/data gmuttiirb/repdbmaker:v1.0 \
  snakemake --configfile config/repdb.yaml --cores <N> --directory /app/data \
  --sdm conda --conda-prefix /conda-envs
```

(`-t` allocates a pseudo-TTY, which is what gets you Snakemake's usual colored
status output — without it everything comes out as plain, uncolored text.)

## Quick start

The workflow is two phases that meet at the **universe** — the enriched table of
every available proteome (id, source, 7 ranks, completeness). Both phases read
the same config file; pick the phase with the target:

| command | produces |
|---------|----------|
| `snakemake sample` | **Pipeline 1 (curation)** — `results/universe/universe.tsv` + `results/universe/repdb.ids` + the taxonomy QC (no sequences fetched) |
| `snakemake build` | **Pipeline 2 (construction)** — the databases: `repdb` and `clusterrepdb` (+ decontamination, stats, `_meta.tsv`) |
| `snakemake` | both (`all`) |

`build` produces the universe first if needed, so it's self-contained; run
`sample` on its own to stop at the universe and review it before building.

To just inspect the available proteomes first:

```bash
snakemake -j 1 --until available_proteomes
# -> results/meta/available_proteomes.tsv, for picking a proteome subset
```

To build everything with the default config:

```bash
snakemake -j 14
```

Producing an actual versioned release — freezing the universe, staging assets,
publishing to GitHub/Zenodo — is a longer process; see [docs/releasing.md](docs/releasing.md).

## Reproducibility

RepDBmaker supports reproducibility at three levels; pick the one your use case
needs.

| Level | What is fixed | How |
|-------|---------------|-----|
| **Parameter** | thresholds, which DBs, which subsets | the config file |
| **Composition** | *which* proteomes and their taxonomy | a pinned **universe** (below) |
| **Artifact** | exact tool builds and, ideally, the exact sequences | conda pin files + Docker digest + a Zenodo deposit of the FASTA |

### Pinned tool versions (conda pin files)

`workflow/envs/*.linux-64.pin.txt` are explicit conda pin files (exact builds
+ md5, `linux-64`) captured from the environments that produced RepDB v1.0. They
pin tools **and** transitive dependencies, and the Docker image is built from
the ones its own envs need (not `benchmark`, see the `Dockerfile`) — but they
aren't Docker-specific: each sits next to its `workflow/envs/<name>.yaml` using
[Snakemake's own pin-file
convention](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#freezing-environments-to-exactly-pinned-packages),
so any `snakemake --sdm conda` run picks it up automatically in place of the
loose `.yaml`, Docker or not. To recreate one directly:

```bash
conda create --prefix ./repdb_homology --file workflow/envs/homology.linux-64.pin.txt
```

The database-building environment pins **DIAMOND 2.1.13**, **MMseqs2 18.8cc5c**,
and **BLAST 2.17.0** (the manuscript versions); the separate benchmark
environment uses **DIAMOND 2.1.21**, which is needed only to read a BLAST
database and is not used to build RepDB.

### Pinned source snapshots

External sources drift, so pin the snapshots under `versions:` in
`config/repdb.yaml`:

- `gtdb` — a specific release (e.g. `release226`), never `latest`.
- `unieuk` — the exported UniEuk taxonomy version.
- `EukProt` — the EukProt version.
- `taxdump` — a dated NCBI taxdump archive (or `latest`).

For sources without stable versioned hosting (UniProt reference-proteome
release, RefSeq virus catalog, P10K), the **universe** below is what actually
freezes them — record the retrieval date alongside your run too.

### Reproducing a release (pinned universe) — two configs

Because several sources are unversioned, re-running the full selection later
yields a *different* set of proteomes. The fix is to freeze the **universe**
and reproduce from it, which splits a release cleanly into two configs:

1. **Curation config** — `config/repdb.yaml`. Generates the universe from the
   live sources (this is the slow, drift-prone part):

   ```bash
   snakemake sample --configfile config/repdb.yaml --sdm conda -j 8
   #   -> results/universe/universe.tsv  +  results/universe/repdb.ids
   ```

2. **Freeze it into a build config** — `make_release` copies the universe (and
   the custom bundle) into `resources/releases/<v>/` and writes
   `resources/releases/<v>/config.yaml`, a self-contained **build config** that
   *pins* those frozen assets:

   ```bash
   snakemake resources/releases/v1/release.yaml
   #   copies universe.tsv (+ custom_bundle.tar.gz) into resources/releases/v1/
   #   and writes resources/releases/v1/config.yaml  with:
   #     dbs: {build: {repdb: {universe: resources/releases/v1/universe.tsv, ...}}}
   ```

3. **Build config** — `resources/releases/<v>/config.yaml`. Rebuild the databases
   from the pinned universe; taxonomy harmonization is skipped, the taxdump is
   the full pinned universe, and the selection re-runs deterministically on it:

   ```bash
   snakemake build --configfile resources/releases/v1/config.yaml --sdm conda -j 8
   ```

Only the sequences are fetched at build time, so this guarantees
composition-level (not bitwise) reproducibility. Commit the light files
(`config.yaml`, `repdb.ids`, `release.yaml`) and publish the heavy assets
(`universe.tsv`, `custom_bundle.tar.gz`) as GitHub-release assets / a Zenodo
deposit; deposit the assembled FASTA too for a bitwise artifact. See
[docs/releasing.md](docs/releasing.md).

`config/example.yaml` is a ready template for the last step, already pointing
`dbs.build.repdb.universe` at `resources/releases/v1/universe.tsv` (download
that asset first) and defining a custom database from your own id list, with
optional clustering / decontamination:

```bash
snakemake build --configfile config/example.yaml --sdm conda -j 8
```

## Troubleshooting

Sometimes things can go wrong while downloading a proteome — rule `db_stats`
fails if any gzipped fasta is malformed, blocking the database fasta until
it's resolved. Full guidance in [docs/troubleshooting.md](docs/troubleshooting.md).

## Documentation

Everything beyond this README lives under [`docs/`](docs/):

- [Choosing an executor](docs/executors.md) — local, SLURM, LSF, site profiles
- [Configuration](docs/configuration.md) — full config reference + eukaryote downsampling
- [Custom databases](docs/custom-databases.md) — adding your own database or proteomes
- [Decontamination](docs/decontamination.md) — the cross-domain contamination filter
- [Outputs and quality control](docs/outputs-and-qc.md) — where everything lands, and the QC reports
- [Utilities and benchmarking](docs/utilities-and-benchmarking.md) — helper scripts, the NR/PhyloDB comparison
- [Troubleshooting](docs/troubleshooting.md)
- [Building & releasing RepDB](docs/releasing.md) — the full release process

## Citation

If you use RepDBmaker or RepDB, please cite:

> Mutti G. and Gabaldón T. Automated reconstruction of reproducible protein
> databases with RepDBmaker. *Protein Science* (2026).
> DOI: <add DOI on acceptance>

## License

See the [LICENSE](LICENSE) file for details.
