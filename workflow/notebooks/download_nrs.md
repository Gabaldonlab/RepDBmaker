# Downloading and Processing NR Databases

This document provides a guide for downloading and processing the NCBI nr and ClusteredNR protein databases.

## Prerequisites

Before starting, ensure you have the following tools installed and available in your PATH:

- `wget` - for downloading files from FTP and HTTP sources
- BLAST tools (`update_blastdb.pl`, `blastdbcmd`) - for managing and querying BLAST databases
- MMseqs2 (`mmseqs`) - for creating and managing sequence databases
- GNU Parallel (`parallel`) - for running commands in parallel
- TaxonKit (`taxonkit`) - for taxonomy manipulation and formatting
- Standard Unix tools (`tar`, `awk`, etc.) - typically pre-installed on Unix-like systems

## NR Database Processing

First, download the taxonomy dump from NCBI

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir taxonomy && tar -xxvf taxdump.tar.gz -C taxonomy
```

Then you can download the nr database using BLAST tools following these steps:

```bash
update_blastdb.pl --num_threads 4 --decompress nr

echo "Verifying md5 checksums..."
ls *.md5 | parallel -j 4 'md5sum -c {} || echo "WARNING: checksum failed for {}"'
ls *.tar.gz | parallel -j 5 'echo "Extracting {}" && tar -xvzf {}'
blastdbcmd -db nr -entry all > nr.fna
blastdbcmd -db nr -entry all -outfmt "%a %T" > nr.fna.taxidmapping
```

Create an MMseqs database from the extracted sequences. Thanks to this MMseqs2 db you can then do taxonomic assignment with nr.

```bash
mmseqs createdb nr.fna nr_mmseqs
mmseqs createtaxdb nr_mmseqs tmp --ncbi-tax-dump taxonomy/ --tax-mapping-file nr.fna.taxidmapping --threads 112
```

Finally, to count the occurrences of each TaxID in the mapping file to analyze taxonomic distribution:

```bash
LC_ALL=C awk -F' ' '{count[$2]++} END {for (i in count) print i"\t"count[i]}' nr.fna.taxidmapping | \
taxonkit reformat -I 1 -f "{C}\t{a}\t{r}\t{d}\t{p}\t{c}\t{o}\t{f}" \
> nr_taxid_counts_tax.tsv
```

## Clustered NR Processing

Downloading ClusteredNR is a bit more tricky and you can follow these steps:

```bash
DEST=./
URL='https://ftp.ncbi.nlm.nih.gov/blast/db/experimental/'

mkdir -p "$DEST"
cd "$DEST"

curl -s "$URL" \
  | grep -oP 'href="\K[^"]+\.tar\.gz' \
  | sort -u \
  | sed "s|^|${URL}|" \
  > to_download.txt

echo "Found $(wc -l < to_download.txt) tarballs."

curl -s "$URL" \
  | grep -oP 'href="\K[^"]+\.md5' \
  | sort -u \
  | sed "s|^|${URL}|" \
  > md5_list.txt

wget -c -i md5_list.txt

cat to_download.txt | parallel -j 4 'wget -c {}'

echo "Verifying md5 checksums..."
ls *.md5 | parallel -j 4 'md5sum -c {} || echo "WARNING: checksum failed for {}"'

ls *.tar.gz | parallel -j 5 'echo "Extracting {}" && tar -xvzf {}'

blastdbcmd -db nr_cluster_seq -entry all > nr_cluster.fna
blastdbcmd -db nr_cluster_seq -entry all -outfmt "%a %T" > nr_cluster.fna.taxidmapping
```

You can then obtain a MMseqs2 taxonomic database following these steps:

```bash
mmseqs createdb nr_cluster.fna nr_cluster_mmseqs
mmseqs createtaxdb nr_cluster_mmseqs tmp --ncbi-tax-dump ../nr/taxonomy/ --tax-mapping-file nr_cluster.fna.taxidmapping --threads 112
```

Finally, in case you are interested, you can count how many sequences are associated for each taxid.

```bash
LC_ALL=C awk -F' ' '{count[$2]++} END {for (i in count) print i"\t"count[i]}' nr_cluster.fna.taxidmapping | \
taxonkit reformat -I 1 -f "{C}\t{a}\t{r}\t{d}\t{p}\t{c}\t{o}\t{f}" \
> clustnr_taxid_counts_tax.tsv
```
