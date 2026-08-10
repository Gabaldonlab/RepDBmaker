# RepDB v1.0 reproducible build image.

FROM condaforge/mambaforge:24.9.2-0@sha256:3b5f55c0b2f447cf94b8f67ca9094f9d2a9d316a4506fd5aa915faaf37ae29d4

# Override at build time: --build-arg SNAKEMAKE_VERSION=...
ARG SNAKEMAKE_VERSION=8.11.6
ARG SLURM_PLUGIN_VERSION=2.7.1
RUN mamba install -y -c conda-forge -c bioconda \
        "snakemake==${SNAKEMAKE_VERSION}" \
        "snakemake-executor-plugin-slurm==${SLURM_PLUGIN_VERSION}"
WORKDIR /app

COPY . /app

RUN snakemake --configfile config/repdb.yaml --sdm conda \
        --conda-prefix /conda-envs --conda-create-envs-only && \
    mamba clean --all -y
