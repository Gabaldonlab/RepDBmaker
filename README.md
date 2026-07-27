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

## Installation

### Requirements

- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (v8.11.6)
- Conda or Miniconda
- Internet access for external downloads
- Sufficient disk space for genome and database files

### Choosing an executor

The workflow ships with several profiles under `workflow/profiles/`. Each cluster
profile names the executor plugin it needs and marks the site-specific values
with `CHANGE_ME`; run any of them with
`snakemake --workflow-profile workflow/profiles/<name> -j <max_parallel_jobs>`.

- `default/` — **local execution, auto-loaded.** Snakemake picks this up
  automatically, so `snakemake -j <cores>` just runs on the current machine with
  no extra flags and no scheduler.
- `slurm/` — **generic SLURM** (`pip install snakemake-executor-plugin-slurm`).
- `lsf/` — **generic LSF** (`pip install snakemake-executor-plugin-lsf`); uses
  the dedicated LSF plugin with native `lsf_queue` / `lsf_project` resources. These are untested profiles so they may require some tailoring.
- `bsc/` — the authors' Barcelona Supercomputing Center profile, kept for
  reference and reproducibility.

The heavy rules request up to 112 cores; on smaller nodes, lower the values in
the profile's `set-threads` block (a commented example is included in each
cluster profile).

If using `--sdm conda`, Snakemake will automatically create the following environments from `workflow/envs/`:

- `workflow/envs/python.yaml` — Python, pandas, matplotlib, polars
- `workflow/envs/homology.yaml` — diamond, mmseqs2, blast
- `workflow/envs/utils.yaml` — taxonkit, csvtk, ncbi-datasets-cli, jq, seqkit
- `workflow/envs/R.yaml` — R and visualization/taxonomy packages

Create all environments without executing the pipeline:

```bash
snakemake --conda-create-envs-only
```

These `.yaml` specs are intentionally loose (minimum bounds only where a feature
requires it) so a fresh install resolves against current packages. To instead
reproduce the **exact** package builds used for the published RepDB v1.0, use the
explicit lockfiles in `workflow/envs/locks/` (see [Reproducibility](#reproducibility)).

### Docker

A Docker image is available at [Docker Hub](https://hub.docker.com/repository/docker/gmuttiirb/repdbmaker/general).

If [Docker](https://docs.docker.com/engine/install/) is installed, run the pipeline from the repository root with:

```bash
docker run --rm -v $(pwd):/app/data gmuttiirb/repdbmaker:v1.0 snakemake --cores 2 --directory /app/data -n
```

For reproducible pulls, reference the image by its immutable digest rather than
the mutable `:v1.0` tag (replace the digest below with the one printed by
`docker inspect --format='{{index .RepoDigests 0}}' gmuttiirb/repdbmaker:v1.0`):

```bash
docker run --rm -v $(pwd):/app/data \
  gmuttiirb/repdbmaker@sha256:5531958bdfe5... \
  snakemake --cores 2 --directory /app/data -n
```

The image builds its conda environments from the pinned lockfiles (see below),
so it reproduces the RepDB v1.0 toolchain rather than re-solving at build time.

## Reproducibility

RepDBmaker supports reproducibility at three levels; pick the one your use case
needs.

| Level | What is fixed | How |
|-------|---------------|-----|
| **Parameter** | thresholds, which DBs, which subsets | the config file |
| **Composition** | *which* proteomes and their taxonomy | a frozen manifest (below) |
| **Artifact** | exact tool builds and, ideally, the exact sequences | conda lockfiles + Docker digest + a Zenodo deposit of the FASTA |

### Pinned tool versions (conda lockfiles)

`workflow/envs/locks/*.linux-64.lock` are explicit conda lockfiles (exact builds
+ md5, `linux-64`) captured from the environments that produced RepDB v1.0. They
pin tools **and** transitive dependencies; the Docker image is built from them.
To recreate one directly:

```bash
conda create --prefix ./repdb_homology --file workflow/envs/locks/homology.linux-64.lock
```

The database-building environment pins **DIAMOND 2.1.13**, **MMseqs2 18.8cc5c**,
and **BLAST 2.17.0** (the manuscript versions); the separate benchmark
environment uses **DIAMOND 2.1.21**, which is needed only to read a BLAST
database and is not used to build RepDB.

### Pinned source databases

External sources drift, so pin the snapshots in `config/repdb.yaml`:

- `gtdb_version` — a specific release (e.g. `release226`), never `latest`.
- `unieuk_version` — the exported UniEuk taxonomy version.

For sources without stable versioned hosting (NCBI taxdump, UniProt release,
RefSeq virus catalog, P10K), record the retrieval date alongside your run.

### Reproducing an exact composition (frozen manifest)

Because several sources are unversioned, running the full selection a year later
yields a *different* set of proteomes. To rebuild the exact composition of a
previous run (e.g. RepDB v1.0) regardless of source drift, freeze its manifest —
`results/taxonomies/repdb_taxonomy.tsv`, a table of each selected ID plus its
seven taxonomic ranks — and feed it back in:

```bash
# 1. after a full run, freeze the selected composition
cp results/taxonomies/repdb_taxonomy.tsv resources/repdb_v1.0.manifest.tsv

# 2. point the config at it
#    dbs:
#      build:
#        repdb:
#          manifest: resources/repdb_v1.0.manifest.tsv

# 3. rebuild — selection and taxonomy harmonization are skipped; only the
#    frozen IDs are fetched and assembled
snakemake -j 14
```

With `manifest` set, the `repdb_taxonomy` selection and the whole taxonomy
harmonization step are bypassed, so the build no longer depends on the live
metadata that decides *which* proteomes are included. The sequences themselves
are still fetched from their sources, so this guarantees composition-level (not
bitwise) reproducibility; deposit the assembled FASTA on Zenodo for a bitwise
artifact.

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

## Building RepDB

The workflow is two phases that meet at the **universe** — the enriched table of
every available proteome (id, source, 7 ranks, completeness). Both phases use the
same `config/repdb.yaml`; you pick the phase with the target:

| command | produces |
|---------|----------|
| `snakemake sample` | **Pipeline 1 (curation)** — `results/universe/universe.tsv` + `results/universe/repdb.ids` + the taxonomy QC (no sequences fetched) |
| `snakemake build` | **Pipeline 2 (construction)** — the databases: `repdb` and `clusterrepdb` (+ decontamination, stats, `_meta.tsv`) |
| `snakemake` | both (`all`) |

`build` produces the universe first if needed, so it is self-contained; run
`sample` on its own to stop at the universe and review it before building.

### First release (v1)

```bash
# 0. put the curated custom proteomes in resources/custom_proteomes/ (CUS<id>.fa|.faa.gz)
snakemake sample --sdm conda -j 8     # universe + selection + QC  -> review
snakemake build  --sdm conda -j 8     # repdb + clusterrepdb
```

`clusterrepdb` is a clustered version of RepDB (the `repdb.ids` selection,
clustered per class) — configured under `dbs.build.custom` in `config/repdb.yaml`.
Since `custom_bundle` there points at a rule output, `build` auto-packages the
release custom bundle from `resources/custom_proteomes/`; unset `custom_bundle`
to read that folder directly instead.

Freeze the run as a reproducible, versioned release:

```bash
snakemake resources/releases/v1/release.yaml   # stages resources/releases/v1/
# then upload the heavy assets and commit the light ones (see resources/releases/README.md)
```

### Reproducing a release / building your own database

Pin a published universe so the harmonization is skipped entirely (fast, drift-
proof) and build whatever subset you want. `config/example.yaml` is a ready
template:

```bash
snakemake build --configfile config/example.yaml --sdm conda -j 8
```

It sets `dbs.build.repdb.universe: resources/releases/v1/universe.tsv` (download
the asset first) and defines a custom db from your own id list, with optional
clustering / decontamination.

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
  gtdb_version: "release226"  # use latest to use most recent version
  unieuk_version: "0.0.1-pre_release"
  type: ["diamond", "mmseqs"]  # Available types: diamond, mmseqs, blastp
  build:
    repdb:
      decontaminate:
        # optimized settings for ContScout benchmarking
        identity: 0.9
        coverage: 0.5
        cov_mode: 3
        prop_euka: 0.5
        filter: soft          # 'hard' removes contaminants; 'soft' flags but keeps them
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

## Eukaryote downsampling

The public eukaryotic sources (UniProt, EukProt, P10K) are heavily redundant, so
before RepDB is built they are thinned by `workflow/scripts/filter_tax.R`. The
three steps are configurable under `dbs.build.repdb.downsample`. The whole block
is optional — omit it, or any individual key, to use the defaults shown below.

```yaml
dbs:
  build:
    repdb:
      downsample:
        # keep only the most complete proteome per species (identical 7-rank
        # lineage). false = keep every proteome.
        remove_duplicated_species: true
        # within each genus, keep at most this many (most complete) species.
        top_n_genuses: 3
        # cap over-represented clades: within each named clade keep at most `n`
        # (most complete) genomes per family; genomes with no resolved family are
        # all kept. `rank` is one of phylum / class / order / family.
        reduce_abundant_clades:
          - {rank: class,  taxon: Opisthokonta, n: 20}
          - {rank: order,  taxon: Ciliophora,   n: 20}
          - {rank: family, taxon: Embryophyta,  n: 20}
```

| Key | Default | Effect |
|-----|---------|--------|
| `remove_duplicated_species` | `true` | Collapse identical species to their single most complete proteome. |
| `top_n_genuses` | `3` | Max number of (most complete) species kept per genus. Genomes with no genus are all kept. |
| `reduce_abundant_clades` | Opisthokonta / Ciliophora / Embryophyta at 20 | Per-family cap applied inside each listed clade. Set to `[]` to disable. |

Notes:

- `reduce_abundant_clades` is a list, so you can add, edit or remove clades
  freely (e.g. `- {rank: phylum, taxon: Chordata, n: 10}`). Within each clade the
  cap is applied **per family**.
- Genomes with no resolved family (an `unassigned_*` family from the taxonomy
  harmonization, or a blank one) are **not** capped — you cannot know where they
  belong taxonomically — so they are all kept.
- Genomes forced in via `clades_to_keep` and all custom (`CUS…`) proteomes bypass
  the downsampling and are always retained.
- A `rank` outside `phylum/class/order/family`, or an entry missing `rank`,
  `taxon` or `n`, aborts the run immediately with a clear error.

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

## Adding custom proteomes

Beyond the public sources, you can inject your own proteomes through the
tab-separated table pointed to by `files.new_genomes` (default
`resources/custom_genomes_repdb.csv`). Each row is **one whole proteome from a
single organism**. The pipeline relies on the following schema, and rows that
violate it are silently dropped or mislabelled rather than rejected — so validate
the table first (see below).

| Column      | Required | Description |
|-------------|----------|-------------|
| `ID`        | yes      | Unique mnemonic, **must start with `CUS`** (e.g. `CUS00001`). Non-`CUS` IDs are silently discarded by `filter_tax.R`. |
| `Species`   | yes      | Organism name (a single organism per row). |
| `Data_type` | yes      | `genome`, `transcriptome`, … (informational). |
| `Fasta`     | yes      | Path to the proteome FASTA for this entry. |
| `Lineage`   | yes      | Exactly **7 `;`-separated ranks** with the prefixes `d__ p__ c__ o__ f__ g__ s__`. |
| `Paper`, `Source`, `Note` | no | Free-text provenance. |

Example row (`Lineage` shown on its own line):

```
CUS00003  Agogonia voluta  transcriptome  /path/CUS00003.fa  <lineage>  <doi>  <url>  <note>
```
```
Lineage = d__Eukaryota;p__Discoba;c__Jakobida;o__Jakobida;f__Ophirinina;g__Agogonia;s__Agogonia voluta
```

Notes on the implicit specification (made explicit here):

- **Single organism, whole proteome.** A row packing proteins from several
  organisms is accepted but all its sequences inherit the row's nominal lineage.
- **Lineages must be internally consistent.** The same taxon name may not appear
  under two different parents across your custom rows (e.g. class `Provora` under
  phylum `Diaphoretickes` in one row and `Metamonada` in another) — such
  contradictions corrupt the taxdump and are rejected.
- **Custom proteomes need not be eukaryotic.** Bacterial/archaeal or viral custom
  proteomes are accepted; each row is checked against the reference matching its
  own domain — Eukaryota vs the harmonized eukaryotic taxonomy, Bacteria/Archaea
  vs GTDB, Viruses vs the RepDB virus taxonomy. The domain (`d__`) must be one of
  `Eukaryota`, `Bacteria`, `Archaea`, `Viruses`.
- **Lineage is checked against the reference for its domain.** A lineage that
  places a *known* taxon under a conflicting parent (e.g. a genus that the source
  taxonomy puts in a different family) is rejected. Genuinely novel taxa (absent
  from the reference) are accepted. A domain with no reference wired in (e.g. in
  reproduce/manifest mode) is validated for schema only.
- **All seven ranks must be present and non-empty** (`d__` through `s__`), because
  the taxdump is built by splitting the lineage into exactly these columns.

### Validating the table

The pipeline validates the table automatically (rule `validate_custom_proteomes`,
gating `custom_taxonomy` and `get_custom_genomes`): a run **aborts** before any
custom lineage is propagated if the table has errors. In a normal run the check
is done against a reference eukaryotic taxonomy built from the public sources
only (`eukaryotes_taxonomy_ref.tsv`, i.e. without the custom proteomes), so the
set of eukaryotes to compare against follows the taxon sampling in your config.
You can also run it standalone:

```bash
Rscript workflow/scripts/check_custom_proteomes.R resources/custom_genomes_repdb.csv
# check lineage conflicts against per-domain reference taxonomies:
Rscript workflow/scripts/check_custom_proteomes.R <table> \
  --reference results/taxonomies/eukaryotes_taxonomy_ref.tsv \
  --reference-prok results/taxonomies/gtdb_taxonomy.tsv \
  --reference-virus results/taxonomies/virus_taxonomy.tsv
# add --check-fasta to also verify each Fasta path exists and is non-empty
```

It **errors** (non-zero exit) on: missing columns, non-`CUS` IDs, duplicate IDs,
lineages that are not exactly seven correctly prefixed non-empty ranks, empty or
duplicate Fasta paths, a domain that is not a recognized superkingdom
(`Eukaryota`/`Bacteria`/`Archaea`/`Viruses`), a `Species` that does not match the
`s__` rank (mislabel / cross-organism row), internal contradictions (a taxon
under conflicting parents across custom rows), and any lineage that conflicts
with the reference taxonomy for its domain (a known taxon placed under a parent it
does not have there). The **only warning** lists the *new coherent lineages* —
custom entries that do not conflict but introduce taxa the reference does not
know — so the novel taxonomy being added can be reviewed before it propagates.

## Decontamination

Adding a `decontaminate:` block under a database enables contamination filtering.
Sequences are clustered across domains; a eukaryotic protein is **flagged** as a
likely contaminant when it is a lone eukaryote in a mixed cluster
(`single_euk`) or when eukaryotes are a minority of its cluster
(`low_euka_prop`, i.e. `euka_prop <= prop_euka`).

```yaml
decontaminate:
  identity: 0.9
  coverage: 0.5
  cov_mode: 3
  prop_euka: 0.5
  filter: soft   # 'soft' (default) or 'hard'
```

- **`filter: soft`** — flagged sequences are **kept** in the database.
- **`filter: hard`** — flagged sequences are **removed** from the database
  (the search indices are built from `<db>_decontaminated.fa.gz`).

**In both modes** a data frame `results/dbs/<db>/decontaminate/contaminants.tsv`
is written, listing every flagged protein with its cluster context: the flag
path, cluster size (`n_clu`), number of eukaryotes (`n_euka`), species count
(`n_species`), `euka_prop`, and the protein's taxonomy. This lets you audit what
was flagged — e.g. whether removals hit HGT candidates or plastid-derived genes —
regardless of the filter mode.

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
- `results/dbs/<db>/decontaminate/contaminants.tsv` — flagged proteins + cluster properties (see [Decontamination](#decontamination))
- `results/dbs/<db>/decontaminate/contaminants.txt` — flagged protein IDs (used by the hard filter)
- `results/dbs/<db>/<db>_decontaminated.fa.gz` — the FASTA the indices are built from (contaminants removed under `filter: hard`)
- `results/dbs/<db>/decontaminate/pair_counts.tsv` — cluster pair counts
- `results/stats/` — database and clustering statistics
- `results/meta/check_resources.txt` — validation of downloaded assets
- `results/taxonomies/` — taxonomy annotations and selection outputs
- `results/taxdump/repdb_taxdump/` — taxonkit taxdump for combined genomes

## Quality control

RepDBmaker emits QC outputs so the heterogeneity of the integrated sources is
visible rather than implicit.

**Taxonomy harmonization report** — rule `taxonomy_harmonization_report` renders
`workflow/notebooks/harmonization_report.Rmd` (produced by default in a normal
run). It compares each proteome's original taxonomy with the UniEuk-harmonized
one and characterizes the cross-framework translation. Outputs to `results/qc/`:

- `harmonization_report.html` — the rendered report: the EukProt backbone tree
  (branches coloured by which database adds them), the summary table, the
  original→UniEuk crosstab heatmaps (NCBI phyla / P10K supergroups → UniEuk
  phyla), and the **ambiguous** proteomes (mapped to more than one UniEuk
  lineage). The notebook also writes these machine-readable tables alongside it:
- `harmonization_uniprot.tsv` / `harmonization_p10k.tsv` — full per-proteome
  before/after detail.
- `harmonization_ambiguous.tsv` — proteomes whose harmonization is contradictory
  (same proteome, two UniEuk lineages), which should be reviewed.

You can also render it standalone from the repo root:

```bash
Rscript -e 'rmarkdown::render("workflow/notebooks/harmonization_report.Rmd")'
```

**Explorer notebook** — `workflow/notebooks/explore_available_proteomes.Rmd` is a
boilerplate for exploring `results/meta/available_proteomes.tsv` (counts by
source/domain/data-type, taxonomic composition, completeness) and for exporting a
selection as an `ids` file for a custom database. Render it from the repo root:

```bash
Rscript -e 'rmarkdown::render("workflow/notebooks/explore_available_proteomes.Rmd")'
```

**Interactive Krona charts** — rule `krona_plot` builds
`results/qc/krona/repdb_krona.html`, a single interactive [Krona](https://github.com/marbl/Krona)
chart with one **selectable dataset per database/taxonomy**: GTDB, NCBI Virus,
EukProt, custom, and — for UniProt and P10K — **both** the source schema and the
UniEuk-harmonized version (`uniprot_ncbi` vs `uniprot_unieuk`, `p10k_native` vs
`p10k_unieuk`). Switching between a source's two datasets shows the harmonization
directly. Open the HTML in a browser.

**Per-organism provenance metadata** — rule `make_db_meta` writes
`results/meta/<db>_meta.tsv` for **every** database it builds (`repdb_meta.tsv`
and one per custom database), one row per organism across **all** sources
(GTDB, NCBI Virus, UniProt, EukProt, P10K, custom). RepDB integrates sources with
different foundations, and the flat FASTA hides that heterogeneity; this table
propagates the intermediate state so it can be filtered on rather than silently
inherited. Columns:

| column | meaning |
|--------|---------|
| `taxid` | the taxid assigned in the taxdump (the exact id the DIAMOND/MMseqs/BLAST indices use) |
| `source_db` | originating database (from superkingdom + id prefix) |
| `taxonomy_authority` | which taxonomy the ranks come from — **GTDB**, **ICTV/NCBI** or **UniEuk**. A `class` (or any rank) is *not* the same concept across these; do not compare ranks across authorities blindly. |
| `k`…`s` | the 7-rank lineage as stored in RepDB |
| `data_type` | genome / transcriptome / single-cell / SAG / MAG … |
| `completeness`, `completeness_metric` | a completeness score **and the metric that produced it** — BUSCO is lineage-specific and CheckM2 is prokaryote-specific, so scores are **only comparable within the same metric** (a 90% BUSCO in a fungus ≠ 90% in a ciliate). Sources with no completeness concept (viruses) are `NA`. |
| `upstream_filter` | the QC/decontamination the **source** applied before RepDBmaker |
| `n_contaminants`, `prop_contaminants` | **RepDBmaker's** own decontamination flagging for the organism (`0` when decontamination was not run) |
| `num_seqs`, `sum_len`, `min_len`, `avg_len`, `max_len` | seqkit size statistics |

Typical uses: restrict completeness comparisons to one `completeness_metric`;
take single-`source_db` subsets when uniform decontamination semantics are needed;
or find organisms from sources whose `upstream_filter` is most aggressive when
looking for over-removed legitimate signal. These make the source heterogeneity
explicit rather than leaving users subject to it.

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

If you use RepDBmaker or RepDB, please cite:

> Mutti G. and Gabaldón T. Automated reconstruction of reproducible protein
> databases with RepDBmaker. *Protein Science* (2026). 
> DOI: <add DOI on acceptance>

## License

See the `LICENSE` file for license details.

## TODOs


* Add snakemake-executor-plugin-slurm to installation instructions
* Test Docker and conda locks
* Test combination of config
* Add QC plots!
* Decide what to do with snakemake version
* Zenodo fasta + taxonomic annotation
* Control check_custom_proteomes.R
* add tests
* better documentation decon
* describe how you fill the missing clades
* test with smk 9
* add package releases
* custom proteomes tar creator (overkill for now)
* change links in tab custom
* custom genomes in available_proteomes

