# RepDB: Comprehensive Protein Sequence Database Builder

RepDB is a Snakemake-based pipeline that builds high-quality, taxonomically-annotated protein sequence databases for homology searches. It combines genomes from multiple public sources (GTDB, UniProt, EukProt, P10K, NCBI RefSeq) and supports custom user-defined databases.

## Features

- 🧬 **Multiple database sources**: GTDB (prokaryotes), UniProt/EukProt/P10K (eukaryotes), NCBI RefSeq (viruses)
- 🔄 **Multiple output formats**: Diamond, MMseqs2 (Blast support planned)
- 📚 **Flexible database selection**: Build RepDB, clustered variants, or fully custom databases
- 🏷️ **Full taxonomic annotation**: Each sequence is mapped to NCBI taxonomy
- 📦 **Scalable**: Easily add custom genomes via configuration
- 🧹 **Quality filtering**: Automatic deduplication and representative selection for eukaryotes

## Quick Start

### 1. Configure
Edit `config/repdb.yaml` to specify which databases to build:
```yaml
dbs:
  type: ["diamond", "mmseqs"]
  build:
    repdb: True              # Standard RepDB
    clustered_repdb: True    # Clustered version
    custom:
      mydb: "/path/to/genomes.tsv"  # Optional custom database
```

### 2. Download Reference Data
```bash
snakemake -j 4 --until online_resources
```

### 3. Build Databases
```bash
snakemake -j 48
```

### Filtering Eukaryotes

The `clades_to_keep.txt` file controls which eukaryotic lineages to include. The pipeline automatically:
- Keeps one genome per genus
- Limits heavily overrepresented families (Opisthokonta, Ciliates) to 20 genomes each
- Removes duplicated species

### Custom Genome Table Format

Custom genome tables must be TSV files with the path to genomes as the first column:
```
/path/to/genome1.faa.gz    d__Eukarya;p__Chordata;c__Mammalia
/path/to/genome2.faa.gz    d__Eukarya;p__Fungi;c__Ascomycota
```

## Running the Pipeline

### Step 1: Download External Resources (internet required)
```bash
snakemake -j 4 --until online_resources
```

### Step 2: Build Databases (compute node)
```bash
snakemake -j 48 -p
```

## Output Structure

```
results/
├── dbs/                          # Built databases
│   ├── repdb/
│   │   ├── repdb.fa.gz          # Full protein sequences
│   │   ├── repdb.map            # Taxonomy mapping
│   │   ├── repdb_diamond/       # Diamond index
│   │   └── repdb_mmseqs/        # MMseqs2 index
│   ├── clustered_repdb/
│   │   ├── repdb.fa.gz          # Clustered sequences
│   │   ├── db_clusters.tsv      # Clustering assignments
│   │   ├── repdb_diamond/
│   │   └── repdb_mmseqs/
│   └── mydb/                     # Custom database (if configured)
│       ├── repdb.fa.gz
│       ├── repdb_diamond/
│       └── repdb_mmseqs/
├── meta/
│   ├── full_genome_table.tsv        # All genomes (unfiltered)
│   ├── repdb_genome_table.tsv       # RepDB genomes (filtered)
│   ├── repdb_meta.tsv               # RepDB statistics
│   ├── clustered_repdb_stats.tsv    # Clustering statistics
│   └── check_resources.txt          # Download status
├── taxonomies/
│   ├── full.tsv                     # All taxonomy annotations
│   ├── repdb.tsv                    # RepDB taxonomy
│   ├── all_eukaryotes.tsv           # All eukaryotic genomes
│   └── eukaryotes.tsv               # Filtered eukaryotes
├── taxdump/
│   └── full_taxdump/                # NCBI taxdump for all genomes
└── log/
    └── dbs/                         # Logs per database

```

## Data Sources

| Source | Type | Coverage |
|--------|------|----------|
| GTDB | Prokaryotes | ~360k bacterial/archaeal genomes |
| UniProt | Eukaryotes | Reference proteomes |
| EukProt | Eukaryotes | High-quality eukaryotic proteomes |
| P10K | Eukaryotes | Plant genomes |
| NCBI RefSeq | Viruses | ~20k viral genomes |
| Custom | Any | User-provided |

## Example Workflows

### Build RepDB + Clustered Version
```yaml
# config/repdb.yaml
dbs:
  type: ["diamond", "mmseqs"]
  build:
    repdb: True
    clustered_repdb: True
```
```bash
snakemake -j 48
```

### Build Multiple Custom Databases
```yaml
# config/repdb.yaml
dbs:
  type: ["diamond", "mmseqs"]
  build:
    repdb: False
    clustered_repdb: False
    custom:
      bacteria: "/data/bacterial_genomes.tsv"
      archaea: "/data/archaeal_genomes.tsv"
      viruses: "/data/viral_genomes.tsv"
```
```bash
snakemake -j 48
# Creates: results/dbs/bacteria/, results/dbs/archaea/, results/dbs/viruses/
```

## Search the Databases

### Using Diamond
```bash
diamond blastp -d results/dbs/repdb/repdb_diamond \
  -q query.fa -o results.txt -f 6 qseqid tseqid pident
```

### Using MMseqs2
```bash
mmseqs search results/dbs/repdb/repdb_mmseqs \
  query_db query_db_results tmp --search-type 2
```

## Disk Space Requirements

| Stage | Space | Notes |
|-------|-------|-------|
| Downloads | ~100GB | Can be deleted after processing |
| Intermediate | ~500GB | GTDB, EukProt FASTA files |
| Final Databases | ~150GB | Diamond + MMseqs2 indices |
| **Total** | **~750GB** | Reusable after cleanup |


## Citation


### Ideas

* Filter gtdb to reduce redundancy? Results were not satisfying.
* genome|gffs db?
