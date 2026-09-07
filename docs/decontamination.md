# Decontamination

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

- **`filter: soft`**: flagged sequences are **kept** in the database.
- **`filter: hard`**: flagged sequences are **removed** from the database
  (the search indices are built from `<db>_decontaminated.fa.gz`).

**In both modes** a data frame `results/dbs/<db>/decontaminate/contaminants.tsv`
is written, listing every flagged protein with its cluster context: the flag
path, cluster size (`n_clu`), number of eukaryotes (`n_euka`), species count
(`n_species`), `euka_prop`, and the protein's taxonomy. This lets you audit what
was flagged (e.g. whether removals hit HGT candidates or plastid-derived genes)
regardless of the filter mode.
