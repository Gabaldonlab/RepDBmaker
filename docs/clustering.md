# Clustering

A `cluster:` block under any database entry (`dbs.build.repdb` or any
`dbs.build.custom.<name>`) is **additive**: it doesn't replace that database,
it makes the pipeline also build a clustered sibling of it, `<db>_clustered`
(e.g. `repdb` + a `cluster:` block gives you both `repdb` and
`repdb_clustered`). Both get their own full set of indices under
`results/dbs/<db>/` and `results/dbs/<db>_clustered/`.

```yaml
dbs:
  build:
    repdb:
      cluster:
        level: class
        identity: 0.9
        coverage: 0.9
        keep_full: True
```

| Key | Default | Effect |
|-----|---------|--------|
| `level` | `class` | The taxonomic rank each clustering run is scoped to. One of `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`. An unrecognized value aborts the run immediately with a clear error. |
| `identity` | `0.9` | Minimum sequence identity for two proteins to land in the same cluster. |
| `coverage` | `0.9` | Minimum alignment coverage for the same. |
| `keep_full` | `true` | `false` tells `cleanup` to treat the unclustered `<db>` as a disposable intermediate (its indices are never built, only `<db>_clustered` ships). `true` keeps both as first-class databases. |

## How it works

Clustering doesn't run once over the whole database. The decontaminated
fasta is first split into one fasta per clade **at the configured `level`**
(the default `class` means: one clade per class present in the database),
using `mmseqs filtertaxseqdb` against the taxdump (or, for genomes with no
resolved taxon at that level, a direct ID lookup instead). Each clade is then
clustered independently with `mmseqs easy-linclust --cluster-mode 2 -e 0.001`
at the configured `identity`/`coverage`, and the representative sequences and
cluster membership from every clade are concatenated back together into the
sibling database's `<db>_clustered.fa.gz` and `<db>_clusters.tsv`.

## Outputs

- `results/dbs/<db>/cluster/cluster_params.yaml`: the resolved clustering
  settings for the run (useful for confirming what a given `<db>_clustered`
  was actually built with).
- `results/dbs/<db>/<db>_clustered.fa.gz`, `results/dbs/<db>/<db>_clusters.tsv`:
  the clustered fasta and cluster membership table, produced as an
  intermediate of `<db>` before being exposed as the sibling database
  `results/dbs/<db>_clustered/`.
