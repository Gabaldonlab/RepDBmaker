# Building & releasing RepDB

How to produce a new release: build the two databases, the **full
decontaminated repdb** and the **clustered repdb**, from a single frozen
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
| **GitHub release** | `universe.tsv`, `custom_bundle.tar.gz` | Pipeline 1 **inputs**, tied to unversioned dbs |
| **Zenodo** | `repdb.fa.gz`, `repdb_clusters.tsv.gz`, `repdb_contaminants.tsv`, `repdb_taxdump.tar.gz`, `repdb_meta.tsv`, stats tables | the data and results, DOI'd and citable |

## Producing a release

Replace `<N>` with the max number of parallel SLURM jobs.

`repdb_clustered` (Pipeline 2's second database) comes from the additive
`cluster:` block already set under `dbs.build.repdb` in `config/repdb.yaml`
(see [Clustering](clustering.md)): it isn't a separate database entry, just
a side effect of that block being present. Put your curated custom proteomes
in `resources/custom_proteomes/` (`CUS<id>.fa|.faa.gz`) before step 1. Since
`custom_bundle` there points at a rule output, `build` auto-packages the
release custom bundle from that folder for you; unset `custom_bundle` to
read the folder directly instead.

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

# 4. STAGE the Zenodo assets - use the PINNED config (same one step 3 built
#    with), not config/repdb.yaml: staging just requests plain paths like
#    results/dbs/repdb/repdb.fa.gz, so Snakemake computes freshness against
#    whatever config you hand it right now. Handed the unpinned config, it
#    reconstructs the DAG as a from-scratch curation run and doesn't trust the
#    already-built outputs, wanting to redo taxonomy harmonization (and,
#    worse, once genome_table.tsv/taxonomy.tsv look freshly regenerated,
#    everything downstream of them too - including the ~86GiB repdb.fa.gz).
#    repdb_clusters.tsv is gzipped during staging (see stage_zenodo_assets):
#    a highly repetitive 2-column TSV that isn't compressed at its source
#    path but shrinks enormously once it is.
snakemake resources/releases/v1/zenodo_stage/zenodo_config.yaml --configfile resources/releases/v1/config.yaml
#    -> resources/releases/v1/zenodo_stage/{repdb.fa.gz, repdb_clusters.tsv.gz,
#       repdb_contaminants.tsv, repdb_taxdump.tar.gz, repdb_meta.tsv,
#       repdb_stats.tsv, repdb_clustered_stats.tsv, SHA256SUMS.txt, zenodo_config.yaml}

# 5. UPLOAD every file in that directory to a new Zenodo deposit. The web UI
#    struggles at this size (repdb.fa.gz alone is tens of GB); use the
#    bucket API instead - get the deposit's `links.bucket` from its API
#    response, then per file:
#      curl --upload-file <file> -H "Authorization: Bearer $ZENODO_TOKEN" \
#        "<bucket-url>/<file>"
#    Also check the draft's "Manage storage" button (Files section) before
#    uploading - the default quota is 50GB per record, and this release
#    needs more; each account has up to 150GB of extra allowance to assign
#    per record without filing a support ticket.
#    Once uploaded, merge zenodo_config.yaml's content (filling in
#    <record_id>) into resources/releases/v1/config.yaml under
#    dbs.build.repdb, and commit that update

# 6. (optional) reclaim disk once both DBs + their indices exist.
#    Set `keep_intermediates: false` in resources/releases/v1/config.yaml, then:
snakemake cleanup --configfile resources/releases/v1/config.yaml
```

## Updating pinned tool versions (conda pin files)

`workflow/envs/<name>.linux-64.pin.txt` pin the exact package builds for each
env. They aren't just what the Docker image is built from: since these are
Snakemake's own native pin files, not bespoke infrastructure (see
[Freezing environments to exactly pinned
packages](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#freezing-environments-to-exactly-pinned-packages)),
any `snakemake --sdm conda` run uses them automatically too, Docker or not.
Snakemake prefers a pin file over its sibling `.yaml` whenever one is present,
falling back to the `.yaml` only if creating from the pin fails.

Regenerate one whenever its `workflow/envs/<name>.yaml` changes, from the real
environment that's actually building the release (not a fresh solve elsewhere:
the point is capturing what actually worked):

```bash
# find the env Snakemake solved for that rule (hash is deterministic - matches
# <name>.yaml's current content; default location shown, adjust if you set a
# custom --conda-prefix):
ls .snakemake/conda/

# dump it:
conda list --explicit --md5 -p .snakemake/conda/<hash> > workflow/envs/<name>.linux-64.pin.txt
```

## Building and publishing the Docker image

```bash
# 1. BUILD - envs are created by `snakemake --conda-create-envs-only` inside
#    the Dockerfile itself (from the pin files above, not a fresh solve; see
#    the Dockerfile's own top-of-file comment for how the addressing works)
docker build -t gmuttiirb/repdbmaker:v1.0 .

# 2. SMOKE-TEST: does the DAG resolve at all?
docker run --rm -t -v $(pwd):/app/data gmuttiirb/repdbmaker:v1.0 \
  snakemake --configfile config/repdb.yaml --cores 2 --directory /app/data \
  --sdm conda --conda-prefix /conda-envs -n

# 3. VERIFY the actual point of this build: baked envs get reused, not
#    recreated. Re-running --conda-create-envs-only should return almost
#    instantly with no "Creating conda environment..." lines - if it prints
#    those, --conda-prefix wasn't passed, or doesn't match /conda-envs, or
#    workflow/envs/*.yaml changed since the image was built.
docker run --rm -t -v $(pwd):/app/data gmuttiirb/repdbmaker:v1.0 \
  snakemake --configfile config/repdb.yaml --directory /app/data \
  --sdm conda --conda-prefix /conda-envs --conda-create-envs-only

# 4. PUSH
docker push gmuttiirb/repdbmaker:v1.0

# 5. GET THE DIGEST and put it in README.md's Docker section
docker inspect --format='{{index .RepoDigests 0}}' gmuttiirb/repdbmaker:v1.0
```

No per-env hash to keep in sync by hand anymore: Snakemake computes and
matches env addresses itself, both at build time (`--conda-create-envs-only`
above) and at run time, from the same `--conda-prefix /conda-envs`, so as
long as every `docker run` also passes `--conda-prefix /conda-envs`, the
addresses always agree by construction. See the Dockerfile's top-of-file
comment for the full mechanism, including the Apptainer/Singularity
trade-off it makes.
