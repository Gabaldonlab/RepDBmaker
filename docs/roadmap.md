# Roadmap and known limitations

Behaviours we know about and have decided **not** to change in the current
release, with what it would take to change them. Each entry records whether the
change would alter the database itself — that is what decides whether it can go
into a patch or has to wait for a new version.

Entries are removed when implemented, and the change noted in the release notes.

---

## 1. `single_euk` counts sequences, not species

**What happens now.** `workflow/scripts/get_contaminants.R` flags a eukaryotic
protein when it is the only eukaryote in a mixed cluster:

```r
single_reps <- cluster_props %>% filter(n_euka == 1) %>% pull(rep)
```

`n_euka` counts *sequences*. As soon as two eukaryotic proteins land in the same
cluster the cluster leaves this path entirely and has to clear
`euka_prop <= prop_euka` instead. Two proteins reach the same cluster either
because the same species is present twice (see item 2) or because one proteome
contributes two paralogs.

**Why it matters.** Detection sensitivity depends on database composition rather
than on the evidence. A cluster of one eukaryotic protein plus one prokaryote is
flagged; the same cluster with a second copy of that eukaryote (2/3 = 0.67) is
not. The effect is one-directional: it can only *hide* contamination, never
invent it, so the reported rates are a floor.

**Possible change.** Count distinct eukaryotic species instead:

```r
n_euka_species = n_distinct(s[k == "Eukaryota"]),
...
single_reps <- cluster_props %>% filter(n_euka_species == 1) %>% pull(rep)
```

Note that `cluster_props` already computes an `n_species`, but it is
`n_distinct(s)` over the *whole* cluster (prokaryotes included) and is never
used — it is not the quantity needed here.

**Cost.** Does **not** change the database under `filter: soft`, because
`decontaminate_db` then symlinks the raw FASTA and the contaminant IDs never
reach the indices. Re-run `get_contaminants` and `make_db_meta` only — minutes.
But the effect is *not* limited to duplicated species: within-proteome paralogs
start being flagged too, so every contamination number moves, probably upward.
Preview before adopting, and expect to update the per-source rates in any text
that quotes them.

After re-running, `snakemake --touch` the usual target: `decontaminate_db`
re-creates its symlink with a fresh mtime and would otherwise make the
clustering and indices look stale.

---

## 2. Force-keeping a clade also exempts it from deduplication

**What happens now.** `filter_tax.R` snapshots every proteome matching
`resources/clades_tokeep.txt` **before** any downsampling, then at the end drops
those mnemos from the downsampled set and re-adds the untouched originals:

```r
manual_keep <- df %>% filter(rowSums(across(everything(), ~. %in% to_keep)) > 0)
...
df <- df %>%
    filter(!mnemo %in% manual_keep$mnemo, !.in_clade) %>%
    bind_rows(reduced, manual_keep, custom_keep)
```

So a force-kept clade bypasses `remove_duplicated_species`, `top_n_genuses` and
`reduce_abundant_clades` alike. `docs/configuration.md` documents the bypass; what
it does not say is that deduplication is included in it.

**Why it matters.** The keep list holds deep-branching, poorly sampled lineages
(CRuMs, Ancyromonadida, Malawimonadida, Telonemia, Rhodelphidia, Anaeramoebidae)
— exactly the clades we want full representation of. But it means duplicate
species survive *only* there, and by item 1 duplicate species suppress the
contamination signal. Retention and detection end up coupled, and the clades the
list protects are the ones whose contamination is most under-reported.

In RepDB v1.0 there are exactly two duplicated eukaryotic species,
*Anaeramoeba ignava* and *A. flamelloides*, both in Anaeramoebidae, each present
once from EukProt and once from UniProt:

```bash
cut -f2-8 results/taxonomies/repdb_taxonomy.tsv | grep Eukaryota | sort | uniq -c | awk '$1>1'
```

**Possible change.** Either apply `remove_duplicated_species` to `manual_keep`
before re-adding it (keeps the clade, drops the redundant copy), or fix item 1 so
that retention no longer affects detection. The second is preferable: keeping
both proteomes of a rare species is defensible, and detection should not care.

**Cost.** Changing what `manual_keep` retains **does** change database
membership: new `repdb.ids`, FASTA, taxdump, rebuilt indices, re-run
decontamination and taxon-aware clustering, and the search benchmark becomes
invalid because it was run against a different database. That is a new version
and new Zenodo artifacts, not a patch.

---

## 3. Per-class contamination summaries are unstable at small n

**What happens now.** Contamination is summarised per eukaryotic class as a
median of per-proteome rates. The highest-ranking classes are represented by one
to four proteomes (Rigifilida and Ancoracysta by one), so the "median" is a
single observation or the average of two, and moves with database composition
rather than with biology.

**Possible change.** Report a pooled rate — total flagged sequences over total
sequences in the class — which is invariant to how the proteomes are split, and
apply a minimum proteome count before ranking classes.

**Cost.** Reporting only, in `workflow/notebooks/eda_report.Rmd`. No compute, no
rebuild.

---

## 4. The release manifest does not record source database versions

**What happens now.** `resources/releases/<version>/config.yaml` snapshots the
`dbs` block, and `release.yaml` records the software version, git commit and
asset checksums. Neither records the `versions:` block — so which GTDB release a
database was built against is not part of the release.

**Why it matters.** `config/default.yaml` currently pins `gtdb: "release226"`, but
the value accepts `latest`, and a release built that way would leave no record of
what `latest` resolved to. Even with a pin, someone reading a released
`config.yaml` cannot see it. This is exactly the ambiguity that made a shifted
genus attribution hard to diagnose: GTDB reshuffles picocyanobacterial genera
between releases, and there was no way to confirm from the release which one had
been used.

**Possible change.** Copy the resolved `versions:` block into
`resources/releases/<version>/config.yaml`, and add the source database versions
to `release.yaml` alongside the software version and git commit. If `latest` was
requested, record what it resolved to, not the literal string `latest`.

**Cost.** Metadata only — a few lines in `make_release`. Does not touch the
database, and can go into the current release.

---

## 5. Contaminant attribution mixes two different units

**What happens now.** `pair_counts.tsv` counts co-occurrences of proteome
**pairs**. Counting its rows per genus — which is what Fig. 3D does — therefore
measures how many *genomes* of that genus share a contaminant cluster with a
proteome. That is a fair measure of how broadly a genus co-occurs, and it is what
the panel is for. But it scales with how many genomes GTDB happens to hold for a
genus, so it cannot be read as a count of contaminant proteins.

**Why it matters.** A single chromatophore-derived protein in *Paulinella
micropora* clustering with homologs from 200 *Cyanobium* and 40
*Prochlorococcus_A* genomes scores 200 versus 40 off that one protein. Any text
quoting "N proteins" from this source is quoting a genome census. Counting
distinct flagged proteins instead is invariant to representation — but it then
turns out to be nearly flat across all Cyanobiaceae genera, because they all
carry homologs of the same genes, which is the honest result: the donor is
identifiable at family level, not genus.

**Possible change.** `eda_prepare.Rmd` already computes both units — genome
counts for the panel (`cooc_by_genus`) and per-protein counts for the two case
studies (`cooc_case`, which reports `Flagged proteins` and `Donor genomes` side
by side). The remaining refinement is to attribute each flagged protein **once**,
to the LCA of the prokaryotic members of its cluster, so the counts partition and
the percentages have an interpretable denominator:

```r
donor_lca <- mixed %>%
  filter(k != "Eukaryota") %>%
  group_by(rep) %>%
  summarise(across(c(k, p, c, o, f, g),
                   ~ if (n_distinct(.x) == 1) .x[1] else NA_character_),
            .groups = "drop")
```

**Cost.** Reporting only, in the notebooks. No rebuild.

---

## 6. Contamination thresholds are knife-edge statistics

**What happens now.** The eukaryotic contamination distribution is dense around
the thresholds used to define "affected" proteomes: 33 proteomes sit between 3%
and 8%, and 60 are above 5%. Fig. 3D selects rows with `contamination >= 5`.

**Why it matters.** Cluster membership shifts slightly between builds — adding
sequences changes which ones linclust picks as representatives — so a handful of
proteomes drift across any fixed line and the count, and the panel height, move
without the biology changing. The count should not be quoted as a result.

**Possible change.** Select figure rows by rank (`slice_max(contamination, n = 40)`)
rather than by threshold. Fixed panel height, comparable across rebuilds, and it
makes the figure about the pattern rather than about a count. Report the source
composition (45 EukProt / 14 P10K / 1 UniProt of the current 60) instead of the
count, since that restates the stable per-source rates.

**Cost.** Reporting only. No rebuild.

---

## 7. Only the BSC profile requests a high-memory node for the decontamination clustering

**What happens now.** `cluster_decontaminate` is the memory bottleneck of a
tree-of-life build: for RepDB v1.0 it needed a MareNostrum 5 GPP-HighMem node
(1,024 GiB). That requirement is declared only in the authors' profile,
`workflow/profiles/bsc/config.yaml`:

```yaml
  cluster_decontaminate:
    runtime: 400
    slurm_extra: "'--qos=gp_bscls' '--constraint=highmem'"
```

The generic `workflow/profiles/slurm/config.yaml` sets `runtime: 400` for the
same rule but no memory request or constraint, and the `lsf` and `default`
profiles say nothing either.

**Why it matters.** A user following `docs/executors.md` with the generic SLURM
profile gets the decontamination clustering scheduled onto whatever a standard
node is on their cluster, where it will be killed at a scale like RepDB's. The
recorded peak resident memory (222 GB, Table S2) does not reveal this either: it
is sampled during execution and sits below a 256 GiB standard node, so the
figure alone suggests the step fits when in practice it does not. Silent
scheduling failures on other clusters are exactly the class of problem worth
closing.

**Possible change.** Declare the requirement portably on the rule itself — a
`resources: mem_mb=...` entry that every executor understands — and keep the
BSC-specific `--constraint=highmem` in the BSC profile as the local translation
of it. Document the bottleneck in `docs/executors.md` and
`docs/decontamination.md` so users on smaller clusters know to expect it before
they hit it.

**Cost.** Scheduling metadata only. Does not change the database.

---

## Pending for the current release

- **`repdb_meta.tsv` needs re-issuing.** Its `n_contaminants` and
  `prop_contaminants` columns are zero for every organism in any copy written
  before the organism-id fix in `analyze_stats.R` (the sequence id is
  `<taxid>_<mnemo>_<accession>`; the tally split on the first field). It ships as
  a Zenodo asset, so a published copy carries the dead columns. Re-run
  `make_db_meta` and replace the asset. `universe.tsv` is unaffected — it has no
  contaminant columns.
