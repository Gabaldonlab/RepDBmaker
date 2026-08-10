# Troubleshooting

Sometimes things can go wrong while downloading a proteome. That is why there is
a step (rule `db_stats`) that will fail if any gzipped fasta is malformed and
will block the creation of the database fasta.

To check for broken files:

```bash
cut -f1 results/dbs/<db>/genome_table.tsv | xargs -I {} sh -c 'gzip -t "{}" || echo "Failed: {}"'
```

Then delete the problematic ones and re-run the pipeline. If the problem
persists, there may be other sorts of problems (the files may be broken or the
current downloading script fails) — we recommend excluding them and finding the
most suitable alternative.
