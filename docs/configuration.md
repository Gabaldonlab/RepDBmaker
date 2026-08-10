# Configuration

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
