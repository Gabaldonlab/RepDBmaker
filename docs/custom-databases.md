# Custom databases

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
  reproduce mode with a pinned universe) is validated for schema only.
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
