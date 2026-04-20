FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="dc3c3317961d4c9fd059d0e6e810e29872fbda7865b3348b98944baf78dfd745"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/R.yaml
#   prefix: /conda-envs/065992a763aefeb76352322080372458
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - bioconductor-biostrings
#     - bioconductor-complexheatmap
#     - bioconductor-ggtree
#     - bioconductor-ggtreeextra
#     - bioconductor-ggmsa
#     - r-adephylo
#     - r-ape
#     - r-circlize
#     - r-data.table
#     - r-ggally
#     - r-ggh4x
#     - r-ggnewscale
#     - r-ggrepel
#     - r-janitor
#     - r-patchwork
#     - r-phytools
#     - r-tidyverse
#     - r-wesanderson
RUN mkdir -p /conda-envs/065992a763aefeb76352322080372458
COPY workflow/envs/R.yaml /conda-envs/065992a763aefeb76352322080372458/environment.yaml

# Conda environment:
#   source: workflow/envs/homology.yaml
#   prefix: /conda-envs/5b69ffa299c6700af906ec1f423eda55
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - diamond
#     - mmseqs2
#     - blast
RUN mkdir -p /conda-envs/5b69ffa299c6700af906ec1f423eda55
COPY workflow/envs/homology.yaml /conda-envs/5b69ffa299c6700af906ec1f423eda55/environment.yaml

# Conda environment:
#   source: workflow/envs/python.yaml
#   prefix: /conda-envs/42bbaac521b8a2ad3736d2332ea6c52b
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - matplotlib
#     - pandas
#     - python
#     - polars
RUN mkdir -p /conda-envs/42bbaac521b8a2ad3736d2332ea6c52b
COPY workflow/envs/python.yaml /conda-envs/42bbaac521b8a2ad3736d2332ea6c52b/environment.yaml

# Conda environment:
#   source: workflow/envs/utils.yaml
#   prefix: /conda-envs/a5d064b61a675a69875701d952d94dbc
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - csvtk
#     - ncbi-datasets-cli
#     - newick_utils
#     - seqkit
#     - taxonkit
#     - jq
#     - wget
#     - tar
#     - unzip
#     - awk
#     - gzip
RUN mkdir -p /conda-envs/a5d064b61a675a69875701d952d94dbc
COPY workflow/envs/utils.yaml /conda-envs/a5d064b61a675a69875701d952d94dbc/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/065992a763aefeb76352322080372458 --file /conda-envs/065992a763aefeb76352322080372458/environment.yaml && \
    mamba env create --prefix /conda-envs/5b69ffa299c6700af906ec1f423eda55 --file /conda-envs/5b69ffa299c6700af906ec1f423eda55/environment.yaml && \
    mamba env create --prefix /conda-envs/42bbaac521b8a2ad3736d2332ea6c52b --file /conda-envs/42bbaac521b8a2ad3736d2332ea6c52b/environment.yaml && \
    mamba env create --prefix /conda-envs/a5d064b61a675a69875701d952d94dbc --file /conda-envs/a5d064b61a675a69875701d952d94dbc/environment.yaml && \
    mamba clean --all -y
