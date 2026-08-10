# Outputs and quality control

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
- `results/dbs/<db>/decontaminate/contaminants.tsv` — flagged proteins + cluster properties (see [Decontamination](decontamination.md))
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
