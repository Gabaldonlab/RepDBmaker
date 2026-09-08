# Getting RepDB v1.0

RepDB v1.0 is deposited on Zenodo under DOI: `<DOI>`. You don't need to
rebuild it from scratch to use it: pick whichever of the two options below
matches what you actually need.

## Just want the data?

If you only need the sequences and their annotations, you can
download the fasta directly from the Zenodo record. No pipeline, no Snakemake, no
setup required.

The record contains:

| File | What it is |
|---|---|
| `repdb.fa.gz` | the protein FASTA |
| `repdb_clusters.tsv.gz` | cluster membership, gzipped (identity + coverage in `docs/clustering.md`) |
| `repdb_contaminants.tsv` | the decontamination report (see `docs/decontamination.md`) |
| `repdb_taxdump.tar.gz` | an NCBI-style taxdump for the sequences in `repdb.fa.gz` |
| `repdb_meta.tsv` | per-organism provenance metadata (see `docs/outputs-and-qc.md`) |
| `repdb_stats.tsv`, `repdb_clustered_stats.tsv` | summary statistics |
| `SHA256SUMS.txt` | checksums for everything above |

After downloading, verify your copy:

```bash
sha256sum -c SHA256SUMS.txt
```

## Want the searchable databases too (DIAMOND / MMseqs2 / BLAST)?

Those aren't hosted on Zenodo: building them locally, from the already-
assembled fasta, is fast enough that shipping multi-hundred-GB index files
isn't worth it. Point RepDBmaker at the record instead, and it fetches the
raw data (checksum-verified) rather than rebuilding it, then builds the
indices on top:

```bash
snakemake build --configfile resources/releases/v1/config.yaml --sdm conda -j <N>
```

`resources/releases/v1/config.yaml` already has this wired in: a pinned
**universe** (so the same 160,045-proteome selection is used, without
re-running taxonomy harmonization) plus the pinned **Zenodo record** (so the
sequences, clusters, contaminants and taxdump are fetched instead of
reassembled). What actually happens when you run this:

1. **Fetched, checksum-verified**: the raw fasta, the cluster table, the
   decontamination report, and the taxdump.
2. **Derived locally, cheaply, no network involved**: the accession maps
   (straight from the fasta headers), the decontaminated fasta (a symlink,
   since RepDB v1.0 uses the soft filter, flagged not removed), and the
   clustered fasta (just the fasta records matching the fetched clusters'
   representative IDs).
3. **Built locally**: the actual DIAMOND, MMseqs2 and BLAST indices, from
   whatever's in your config's `dbs.type`.

Only step 3 does any real computation. Steps 1 and 2 are what make this
much faster than a full rebuild.

Docker works the same way; see [Docker](../README.md#docker) for the flags
it needs, and swap in `resources/releases/v1/config.yaml` as the
`--configfile`.

## Which one reproduces exactly what we published?

Both do, in the sense that matters: the sequences, clusters, and
decontamination report you get are byte-identical to the ones behind the
manuscript's figures and statistics either way, since they're the same
files, checksummed. The only difference is whether you also build search
indices locally. See [Reproducibility](../README.md#reproducibility) in the
main README for how this compares to the other reproducibility levels
RepDBmaker supports (rebuilding from a pinned universe, or from live
sources entirely).
