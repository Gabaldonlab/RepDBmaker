# Choosing an executor

The workflow ships with several profiles under `workflow/profiles/`. Each cluster
profile names the executor plugin it needs and marks the site-specific values
with `CHANGE_ME`; run any of them with
`snakemake --workflow-profile workflow/profiles/<name> -j <max_parallel_jobs>`.

- `default/`: **local execution, auto-loaded.** Snakemake picks this up
  automatically, so `snakemake -j <cores>` just runs on the current machine with
  no extra flags and no scheduler.
- `slurm/`: **generic SLURM** (`pip install snakemake-executor-plugin-slurm`).
- `lsf/`: **generic LSF** (`pip install snakemake-executor-plugin-lsf`); uses
  the dedicated LSF plugin with native `lsf_queue` / `lsf_project` resources.
  These are untested profiles so they may require some tailoring.
- `bsc/`: the authors' Barcelona Supercomputing Center profile, kept for
  reference and reproducibility.

The heavy rules request up to 112 cores; on smaller nodes, lower the values in
the profile's `set-threads` block (a commented example is included in each
cluster profile).
