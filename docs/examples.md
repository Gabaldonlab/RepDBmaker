# Example commands

A cookbook of runnable commands for common tasks, grouped by what you're
trying to do. All assume the repository root as the working directory; add
`--sdm conda` (and, inside Docker, `--conda-prefix /conda-envs`) to any of
these to have Snakemake manage the conda envs itself, as in
[Installation](../README.md#installation).


## Try the whole pipeline fast, on a tiny subset

`test: true` (or an integer) in the config keeps only a handful of genomes
per source, so a full `sample` + `build` run finishes in minutes instead of
hours.

```bash
snakemake --configfile config/repdb.yaml --config test=true --sdm conda -j 4
```

The large downloads (GTDB tarball, EukProt archive, virus zip) and the QC/Krona
reports still run on the full sources even in test mode.

## Inspect available proteomes before selecting anything

```bash
snakemake --configfile config/repdb.yaml -j 1 --until available_proteomes
```

Useful for browsing what's available (counts by source/domain/completeness)
before deciding on a subset; see the explorer notebook in
[Outputs and quality control](outputs-and-qc.md).

## Build everything with the default config

```bash
snakemake --configfile config/repdb.yaml --sdm conda -j <N>
```

Runs both phases (`sample` then `build`) against `config/repdb.yaml`. Use
`snakemake sample --configfile config/repdb.yaml --sdm conda -j <N>` alone
to stop at the universe and review the selection first.

## Build with a different config file

```bash
snakemake build --configfile path/to/custom.yaml --sdm conda -j <N>
```

See [Configuration](configuration.md) for the full config reference.

## Build only a custom database, no RepDB

`snakemake build` builds every database your config defines: `repdb` if
`dbs.build.repdb` is set, plus every entry under `dbs.build.custom`. Point at
a config that defines only a custom entry (no `repdb:` block) to build just
that:

```bash
snakemake build --configfile config/example.yaml --sdm conda -j <N>
```

See [Custom databases](custom-databases.md) for the config schema.

## Run on a scheduler instead of locally

```bash
snakemake --configfile config/repdb.yaml --workflow-profile workflow/profiles/slurm -j <N>
```

`<N>` is the max number of jobs submitted in parallel, not cores per job (set
per rule in the profile). See [Choosing an executor](executors.md) for the
other available profiles.

## Reproduce a published release exactly (composition-level)

```bash
snakemake build --configfile resources/releases/v1/config.yaml --sdm conda -j <N>
```

Skips taxonomy harmonization entirely and re-derives the selection
deterministically from a frozen universe. See
[Reproducing a release](../README.md#reproducibility) and
[Building & releasing RepDB](releasing.md).

## Validate a custom proteomes table before running anything

```bash
Rscript workflow/scripts/check_custom_proteomes.R resources/custom_genomes_repdb.csv
```

Catches schema and lineage-conflict errors up front rather than mid-run. See
[Custom databases](custom-databases.md#validating-the-table) for the full
flag set (checking against a specific reference, verifying FASTA paths exist).

## Extract a sequence subset from a built database

```bash
python workflow/scripts/get_fasta.py <id_file> <output_fasta> results/dbs/<db>/<db>_mmseqs
```

See [Utilities and benchmarking](utilities-and-benchmarking.md).

## Render a QC report standalone

```bash
Rscript -e 'rmarkdown::render("workflow/notebooks/harmonization_report.Rmd")'
```

Useful for re-rendering after tweaking the notebook without re-running the
whole pipeline. See [Outputs and quality control](outputs-and-qc.md) for what
it covers and for the explorer notebook's equivalent command.

## Reclaim disk once a build is finished

Set `keep_intermediates: false` in your config, then:

```bash
snakemake cleanup --configfile <your-config>.yaml
```

Deletes the big non-final intermediates (raw and full-decontaminated fastas,
clustering scratch, decontamination scratch) once every database's final
indices exist. Never runs as part of `all`/`build`, and is a safe no-op if
`keep_intermediates` is left at its default (`true`).


## Run the NR/PhyloDB annotation benchmark

Needs its own config, `config/benchmark.yaml`, edited first to point at
reference databases and a query set on your own filesystem (paths in the
shipped file are institution-specific placeholders):

```bash
snakemake -s workflow/rules/benchmark.smk --configfile config/benchmark.yaml --sdm conda -j <N>
```

See [Utilities and benchmarking](utilities-and-benchmarking.md).
