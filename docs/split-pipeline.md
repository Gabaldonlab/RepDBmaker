# Splitting RepDBmaker into two pipelines

Status: **in progress** (branch `split-pipe`). Stage 1 (the seam) is implemented;
Stages 2–3 (physically moving the build rules into their own module and pivoting
the taxdump/fetch onto the universe) are planned below.

## Why

RepDBmaker does two very different jobs in one workflow:

1. **Curation / sampling** — download source *metadata*, harmonize taxonomy
   (UniEuk/EukProt — the fragile, sometimes-unretrievable step), and choose an
   optimal taxon sampling across the tree of life. Metadata-light, fragile,
   changes ~once per release.
2. **Construction** — fetch the *sequences* for a chosen set of ids and build the
   FASTA, taxdump and search indices, with clustering/decontamination. Sequence-
   heavy, mechanical, rebuilt whenever you want a new variant.

Keeping them together means every build re-touches the fragile harmonization, the
reproduce path has to be a special `manifest` mode, and RepDB vs custom databases
are two code paths for what is really one operation ("build a db from a set of
ids"). Splitting removes all three problems.

## The seam: the `universe`

The two pipelines meet at one enriched, versionable artifact:

`results/universe/universe.tsv` — one row per **available** proteome (every
source), with the columns Pipeline 2 needs so it never re-harmonizes:

| column | meaning |
|--------|---------|
| `mnemo` | organism id (also the fetch key) |
| `source_db` | gtdb / virus / uniprot / eukprot / p10k / custom |
| `k`…`s` | the 7 harmonized ranks (feed the taxdump and the selection) |
| `data_type` | genome / transcriptome / single-cell / SAG / MAG … |
| `completeness` | BUSCO/CheckM2 where one exists (euks + gtdb), else `NA` |
| `annotated` | `FALSE` for known-unannotated proteomes (dropped by selection) |

Plus `results/universe/repdb.ids` — the selected RepDB composition (all
prokaryotes + all viruses + the selected eukaryotes).

Freezing `universe.tsv` per release freezes **taxonomy + BUSCO**, so Pipeline 2
reproduces a release deterministically (same taxdump — taxids are hashed
per-lineage by `taxonkit create-taxdump`, so they're set-independent — and, with
the same code+config, the same selection). The universe also doubles as the
"available proteomes" menu and as the base for the per-db provenance table.

```
   Pipeline 1 (sampling/curation)              Pipeline 2 (construction)
   ────────────────────────────               ─────────────────────────
   metadata downloads                         universe.tsv  ─┐
   harmonization (get_tax.R)      universe    repdb.ids     ─┤→ taxdump (full)
   filter_tax (downsampling)   ═══════════▶   custom ids     │→ fetch sequences
   build_universe  ────────────▶ universe.tsv                │→ FASTA → indices
   emit_repdb_ids  ────────────▶ repdb.ids                   │→ cluster/decon/meta
```

## Design decisions (agreed)

- **Single-file** universe (ranks + completeness together), not two files.
- **Keep `repdb.ids`** and **report** how a fresh selection differs from it
  (`rule repdb_ids_diff` → `results/qc/repdb_ids_diff.txt`) — a **soft,
  non-fatal** provenance diff, never a build gate. Configure the reference with
  `files.repdb_ids_reference`.
- **Tie-break determinism**: `filter_tax` now `arrange(mnemo)`s before the
  `slice_max` steps so equal-BUSCO ties don't depend on input file order.
- **Custom proteomes live in Pipeline 2**: a build extends the effective universe
  with validated user additions (`check_custom_proteomes` already validates each
  row against the reference for its domain). Pipeline 1 need not know about them.
- **Modules in one repo** via Snakemake `module`s, shared `scripts/` + `envs/`,
  one config with `sampling:` and `build:` sections, and a top-level target that
  runs both end-to-end for the simple user.

## Stage 1 — DONE (this branch)

- `workflow/scripts/build_universe.R` + `rule build_universe` → `universe.tsv`.
- `rule emit_repdb_ids` → `repdb.ids`; `rule repdb_ids_diff` → soft diff.
- `filter_tax` tie-break (`arrange(mnemo)`).
- Wired into `rule all` (normal mode); unit-tested (`.tests/unit/test_build_universe.R`).

These are additive — the existing pipeline is unchanged and still runs.

## Stage 2 — split the rules into two modules

Move rules into `workflow/rules/sample/` and `workflow/rules/build/`, expose each
as a Snakemake `module`, and connect them by the universe file paths (not
`rules.X` cross-refs). Assignment:

- **sample:** `download_*`/`get_*` **metadata** rules, `*_taxonomy`
  (harmonization), `eukaryotes_taxonomy_ref`, `select_repdb_eukaryotes`,
  `build_universe`, `emit_repdb_ids`, `krona_*`, `taxonomy_harmonization_report`.
- **build:** sequence downloads (gtdb tarball, eukprot tgz, virus zip), the
  per-proteome `get_*_genomes`/`extract_virus`, `create_full_taxdump` (→ taxdump
  from the universe), `make_db_*`, `cluster.smk`, `decontaminate.smk`,
  `db_stats`, `make_db_meta`, custom validation + `custom_dbs_taxonomy`.

Two couplings to resolve during the move (both currently metadata-derived-from-
sequences or live-metadata):

1. `gtdb_species_clusters` derives `gtdb.ids` by *extracting the tarball*. In the
   split, the gtdb ids come from the universe (the gtdb rows), and the tarball is
   only extracted for sequences.
2. Fetch rules that read live source metadata (e.g. the UniProt taxid lookup)
   must derive everything from `id + archive`; anything else moves into the
   universe or is dropped.

## Stage 3 — pivot reproduce mode onto the universe

Replace the `manifest` special-casing: Pipeline 2 takes a `universe:` input that
is either generated by Pipeline 1 or a pinned `resources/repdb_vX.universe.tsv`.
`create_full_taxdump` builds the **full** taxdump from the universe in both cases
(fixing the manifest taxdump-shrink), and the selection is re-run (deterministic)
or read from a pinned `repdb.ids`. `manifest`/`REPDB_MANIFEST` then retire.
