# Building & releasing RepDB (v1)

How to build the two v1 databases — the **full decontaminated repdb** and the
**clustered repdb** — reproducibly from a single frozen universe.

## What v1 contains

A single `cluster` block on `repdb` makes the pipeline emit **two first-class
databases** from one entry:

| DB | path | what it is |
|----|------|------------|
| `repdb` | `results/dbs/repdb/` | full RepDB, **decontaminated (soft)** + diamond/mmseqs indices |
| `repdb_clustered` | `results/dbs/repdb_clustered/` | repdb's decontaminated fasta **clustered per class** + its own indices |

Both are built from the **same** frozen universe, so the paper figures and the
released databases describe the same data.

## Clustering is an additive variant

`cluster:` on a db does **not** replace its final with the clustered version — it
adds a parallel sibling db `<db>_clustered`. Key properties:

- The clustered variant is built by **clustering the parent's final fasta**, i.e.
  the *decontaminated* repdb (chain: `raw → decontaminate → cluster`). So it
  **inherits decontamination automatically** — there is no separate, error-prone
  decontamination of the clustered set, and no `ids:` coupling to keep in sync.
- Because it derives from repdb's own output, it also **inherits test mode**: a
  `test: N` run subsamples repdb, and the clustered variant follows.
- The same `cluster` block works on any **custom** db.

With `filter: soft` (the default), decontamination removes no sequences — it only
writes `contaminants.tsv`. So the clustered variant's sequences *are* the
decontaminated sequences, and the authoritative contaminant report is repdb's
`results/dbs/repdb/decontaminate/contaminants.tsv` (the clustered representatives
keep the same `taxid_code_AAxxx` IDs, so it maps straight onto them).

> If a future release should **physically remove** contaminants (`filter: hard`),
> nothing about the wiring changes — the clustered variant already clusters the
> hard-decontaminated fasta, because it consumes the parent's final output.

## Commands

One frozen universe → freeze v1 → build both DBs. Replace `<N>` with the max
number of parallel SLURM jobs.

```bash
# 0. In config/repdb.yaml make sure test mode is OFF:
#      test: false
#    (a release is never frozen in test mode; make_release also forces it off)

# 1. CURATION — build the frozen universe + RepDB selection ONCE
snakemake sample --configfile config/repdb.yaml \
  --workflow-profile workflow/profiles/bsc -j <N>
#    -> results/universe/universe.tsv
#    -> results/universe/repdb.ids

# 2. FREEZE v1 — pin the universe (+ custom bundle) into a self-contained config
snakemake resources/releases/v1/release.yaml --configfile config/repdb.yaml
#    -> resources/releases/v1/config.yaml   (pins universe; carries the cluster block)
#    -> resources/releases/v1/universe.tsv, repdb.ids, release.yaml, custom_bundle.tar.gz
#    prints the `gh release create` + `git add/commit` commands for publishing

# 3. BUILD both DBs from the pinned config (reproduce mode: skips harmonization,
#    re-runs the selection deterministically on the frozen universe)
snakemake build --configfile resources/releases/v1/config.yaml \
  --workflow-profile workflow/profiles/bsc -j <N>
#    -> results/dbs/repdb/repdb_{diamond,mmseqs}
#    -> results/dbs/repdb_clustered/repdb_clustered_{diamond,mmseqs}

# 4. (optional) reclaim disk once both DBs + their indices exist.
#    Set `keep_intermediates: false` in resources/releases/v1/config.yaml, then:
snakemake cleanup --configfile resources/releases/v1/config.yaml
```

The relevant `config/repdb.yaml` block:

```yaml
dbs:
  build:
    repdb:
      decontaminate:
        filter: soft
        # ... benchmarked params ...
      cluster:            # additive: also emits repdb_clustered
        level: class
        identity: 0.9
        coverage: 0.9
        keep_full: true   # false -> cleanup drops the full repdb; only repdb_clustered ships
```

## Notes

- **One universe for draft + release.** The draft numbers must come from this
  same build (`resources/releases/v1/`). Do not freeze a second universe, or the
  draft would describe a database you did not release.
- **Rebuilding v1 later** is `snakemake build --configfile
  resources/releases/v1/config.yaml ...` again — the pinned universe reproduces
  the exact composition regardless of how the online sources have drifted.
- **`keep_full`** (per db, under `cluster:`): `true` (default) ships both the full
  and clustered dbs; `false` treats the full db as a disposable intermediate — its
  indices are never built and `cleanup` removes its fasta, leaving only
  `<db>_clustered`.
- **`cleanup` is a deliberate double opt-in**: never part of `all`/`build`, and a
  no-op unless `keep_intermediates: false`. It deletes only regenerable
  intermediates (raw + full-decontaminated fastas, the clustering mmseqs DB and
  per-clade scratch, decon detection scratch, and the parent-side copy of the
  clustered fasta that was hard-linked into the sibling db). It keeps each DB's
  final fasta, indices, maps and provenance TSVs. It also removes the raw
  `repdb.fa.gz` (the decon reference for custom dbs), so run it only after **all**
  DBs are built.

## Fast path: reproducing from a Zenodo-pinned release

Once v1's heavy assets are on Zenodo, `snakemake build` can fetch them instead
of re-deriving repdb from raw sources - a second, further-downstream seam than
step 3 above. Step 3 still redoes assembly, decontamination-detection and
taxon-aware clustering deterministically on the pinned universe (hours); this
fetches their outputs directly (mostly bandwidth-bound), and only
`diamond makedb`/`mmseqs createdb` still run locally, so the indices always
match whatever diamond/mmseqs version is actually installed rather than
shipping a version-pinned binary index.

```bash
# after step 3 has built repdb once, stage + checksum its heavy assets and
# print a ready-to-paste `zenodo:` config block:
snakemake resources/releases/v1/zenodo_stage/MANIFEST.txt --configfile config/repdb.yaml

# upload everything in resources/releases/v1/zenodo_stage/ to a new Zenodo
# deposit (web UI or the REST API - no official CLI), then paste the printed
# block into resources/releases/v1/config.yaml under dbs.build.repdb:
#   zenodo:
#     base_url: "https://zenodo.org/records/<record_id>/files"
#     files:
#       fasta: {name: repdb.fa.gz, sha256: "..."}
#       ...

# reproduce repdb + repdb_clustered from Zenodo instead of rebuilding:
snakemake results/dbs/repdb/repdb_diamond results/dbs/repdb/repdb_mmseqs \
  results/dbs/repdb_clustered/repdb_clustered_diamond \
  results/dbs/repdb_clustered/repdb_clustered_mmseqs \
  --configfile resources/releases/v1/config.yaml -j <N>
```

Notes:

- **Only one fasta is hosted.** `repdb_clustered.fa.gz` is not a separate
  Zenodo asset - it's derived locally from the raw fasta + `repdb_clusters.tsv`
  by extracting the cluster tsv's representative IDs (column 1, mmseqs' own
  convention) with `seqkit grep`, instead of re-running `mmseqs easy-linclust`.
  Both databases reproduce from one download.
- **`filter: soft` only** (v1.0's default - flags contaminants, removes
  nothing). The decontaminated fasta is then a symlink of the raw one, same as
  `decontaminate_db`'s own soft-mode behaviour. A `filter: hard` release would
  need its already-filtered fasta hosted under its own key instead.
- **Requesting `build` directly still touches Pipeline 1.** The QC/provenance
  tables (`results/stats/*_stats.tsv`, `results/meta/repdb_meta.tsv`) are
  genome_table/universe-derived and aren't part of this fast path - only the
  four index targets above are. `snakemake sample` isn't needed for those.
- Every fetch is checksum-verified against the sha256 pinned in the config
  (not Zenodo's API at run time), so the release config is self-verifying.
  See `workflow/rules/zenodo_fetch.smk` and `repdb_zenodo()` in `taxonomy.smk`.
