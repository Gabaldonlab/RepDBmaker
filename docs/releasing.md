# Building & releasing RepDB

How to produce a new release: build the two databases — the **full
decontaminated repdb** and the **clustered repdb** — from a single frozen
universe, and publish everything. Replace `v1` below with the new version.

## What a release contains

| DB | path | what it is |
|----|------|------------|
| `repdb` | `results/dbs/repdb/` | full RepDB, **decontaminated (soft)** + diamond/mmseqs indices |
| `repdb_clustered` | `results/dbs/repdb_clustered/` | repdb's decontaminated fasta **clustered per class** + its own indices |

Both are built from the **same** frozen universe, so the paper figures and the
released databases describe the same data.

## Where things live

Three hosting tiers, split by what each artifact *is* rather than by size:

| tier | artifacts | why |
|---|---|---|
| **git** (this repo) | `config.yaml`, `repdb.ids`, `release.yaml` | versioned code |
| **GitHub release** | `universe.tsv`, `custom_bundle.tar.gz` | Pipeline 1 **inputs** - tied to unversioned dbs |
| **Zenodo** | `repdb.fa.gz`, `repdb_clusters.tsv`, `repdb_contaminants.tsv`, `repdb_taxdump.tar.gz`, `repdb_meta.tsv`, stats tables | the data and results - DOI'd and citable |

## Producing a release

Replace `<N>` with the max number of parallel SLURM jobs.

```bash
# 0. config/repdb.yaml: test mode off (test: false). A release is never
#    frozen in test mode; make_release also forces it off regardless.

# 1. CURATION - build the frozen universe + RepDB selection ONCE
snakemake sample --configfile config/repdb.yaml \
  --workflow-profile workflow/profiles/bsc -j <N>

# 2. FREEZE - pin the universe (+ custom bundle) into a self-contained config
snakemake resources/releases/v1/release.yaml --configfile config/repdb.yaml
#    -> resources/releases/v1/{config.yaml, repdb.ids, release.yaml,
#       universe.tsv, custom_bundle.tar.gz}
#    prints the `gh release create` + `git add/commit` commands - run them now

# 3. BUILD both DBs from the pinned config (reproduce mode: skips
#    harmonization, re-runs the selection deterministically on the frozen universe)
snakemake build --configfile resources/releases/v1/config.yaml \
  --workflow-profile workflow/profiles/bsc -j <N>

# 4. STAGE the Zenodo assets
snakemake resources/releases/v1/zenodo_stage/zenodo_config.yaml --configfile config/repdb.yaml
#    -> resources/releases/v1/zenodo_stage/{repdb.fa.gz, repdb_clusters.tsv,
#       repdb_contaminants.tsv, repdb_taxdump.tar.gz, repdb_meta.tsv,
#       repdb_stats.tsv, repdb_clustered_stats.tsv, SHA256SUMS.txt, zenodo_config.yaml}

# 5. UPLOAD every file in that directory to a new Zenodo deposit (web UI or
#    the REST API - no official CLI), then merge zenodo_config.yaml's content
#    (filling in <record_id>) into resources/releases/v1/config.yaml under
#    dbs.build.repdb, and commit that update

# 6. (optional) reclaim disk once both DBs + their indices exist.
#    Set `keep_intermediates: false` in resources/releases/v1/config.yaml, then:
snakemake cleanup --configfile resources/releases/v1/config.yaml
```
