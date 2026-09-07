# Utilities and benchmarking

## Utilities

The repository includes helper scripts for working with generated databases.

- `workflow/scripts/get_fasta.py`: extract a subset of FASTA sequences from an MMseqs2 database

Example usage:

```bash
python workflow/scripts/get_fasta.py <id_file> <output_fasta> <db_mmseqs>
```

## Benchmarking

A benchmark workflow is available in `workflow/rules/benchmark.smk` and uses
`config/benchmark.yaml` to compare RepDB against reference databases such as NR
and clustered NR. It's a separate entry point (its own Snakefile/config, not
part of building RepDB itself), with its own conda env
(`workflow/envs/benchmark.yaml`) pinning **DIAMOND 2.1.21**, a newer build than
the one RepDB is built with (2.1.13), needed only to read a BLAST database and
not used anywhere in building RepDB. See the NB in the `Dockerfile` for how to
bake this env if you need it.

A results notebook is available at `workflow/notebooks/comparison.Rmd`.
