# RepDB v1.0

A Tree-of-Life-scale protein sequence database: about 482 million protein
sequences from 160,045 genomes and transcriptomes, spanning prokaryotes
(GTDB), viruses (NCBI Virus), and a taxonomically balanced selection of
eukaryotes (EukProt, the Protist 10000 Genomes Project, UniProt Reference
Proteomes, and a curated custom set). Taxonomy is standardized across
sources: GTDB for prokaryotes, ICTV/NCBI for viruses, and the UniEuk
framework for eukaryotes.

Built with [RepDBmaker](https://github.com/Gabaldonlab/RepDBmaker), a
Snakemake pipeline for constructing reproducible protein databases like
this one.

## What's in this record

| File | What it is |
|---|---|
| `repdb.fa.gz` | the protein FASTA, headers formatted `<taxid>_<mnemonic>_<accession>` |
| `repdb_clusters.tsv.gz` | cluster membership, gzipped (class-level, 90% identity/coverage) |
| `repdb_contaminants.tsv` | the decontamination report, every flagged sequence with its cluster context |
| `repdb_taxdump.tar.gz` | an NCBI-style taxdump covering every taxon in `repdb.fa.gz` |
| `repdb_meta.tsv` | per-organism metadata: source database, taxonomy, completeness, contamination stats |
| `repdb_stats.tsv`, `repdb_clustered_stats.tsv` | sequence-count and length summaries |
| `SHA256SUMS.txt` | checksums for every file above |

RepDB v1.0 uses the *soft* decontamination filter: flagged sequences are
kept in `repdb.fa.gz`, not removed, and marked in `repdb_contaminants.tsv`
instead.

**Not included, because it can be derived from what is**: the
clustered fasta itself. `repdb_clusters.tsv.gz`'s first column is the
representative sequence ID for every cluster, so it's exactly the set of
headers to keep:

```bash
zcat repdb_clusters.tsv.gz | cut -f1 | sort -u > cluster_reps.txt
seqkit grep -f cluster_reps.txt repdb.fa.gz -o repdb_clustered.fa.gz
```

## First, verify your download

```bash
sha256sum -c SHA256SUMS.txt
```

## Using this data directly, no pipeline needed

You don't need RepDBmaker or Snakemake to build search databases from
these files. First, extract the taxdump and derive an accession-to-taxid
map from the fasta headers (both diamond and mmseqs need one; BLAST wants
one without the header line):

```bash
mkdir repdb_taxdump && tar -xzf repdb_taxdump.tar.gz -C repdb_taxdump

echo -e "accession.version\ttaxid" > repdb.map
zcat repdb.fa.gz | grep '^>' | sed 's/^>//' | awk -F'_' '{print $0"\t"$1}' >> repdb.map
awk 'NR>1' repdb.map > repdb_nohead.map
```

**DIAMOND:**
```bash
diamond makedb --in repdb.fa.gz -d repdb_diamond \
  --taxonnodes repdb_taxdump/nodes.dmp --taxonmap repdb.map
```

**MMseqs2:**
```bash
mmseqs createdb repdb.fa.gz repdb_mmseqs
mmseqs createtaxdb repdb_mmseqs tmp \
  --ncbi-tax-dump repdb_taxdump --tax-mapping-file repdb_nohead.map
```

**BLAST:**
```bash
gunzip -c repdb.fa.gz | makeblastdb -in - -parse_seqids -taxid_map repdb_nohead.map \
  -dbtype prot -out repdb_blastp -title repdb
```

## Using it through RepDBmaker instead

If you'd rather have this fetched, verified, and indexed automatically,
point [RepDBmaker](https://github.com/Gabaldonlab/RepDBmaker) at this
record: it downloads exactly the four files above, checksums them, derives
the clustered fasta and the maps the same way shown here, and builds all
three index types in one command. See
[Getting RepDB v1.0](https://github.com/Gabaldonlab/RepDBmaker/blob/main/docs/get-repdb.md)
in the repository for the exact steps.
