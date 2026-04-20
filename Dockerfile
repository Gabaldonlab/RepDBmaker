FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="dc3c3317961d4c9fd059d0e6e810e29872fbda7865b3348b98944baf78dfd745"


RUN mamba install -y -c conda-forge -c bioconda snakemake snakemake-executor-plugin-slurm
WORKDIR /app

RUN mkdir -p /conda-envs/065992a763aefeb76352322080372458
COPY workflow/envs/R.yaml /conda-envs/065992a763aefeb76352322080372458/environment.yaml

RUN mkdir -p /conda-envs/5b69ffa299c6700af906ec1f423eda55
COPY workflow/envs/homology.yaml /conda-envs/5b69ffa299c6700af906ec1f423eda55/environment.yaml

RUN mkdir -p /conda-envs/42bbaac521b8a2ad3736d2332ea6c52b
COPY workflow/envs/python.yaml /conda-envs/42bbaac521b8a2ad3736d2332ea6c52b/environment.yaml

RUN mkdir -p /conda-envs/a5d064b61a675a69875701d952d94dbc
COPY workflow/envs/utils.yaml /conda-envs/a5d064b61a675a69875701d952d94dbc/environment.yaml

# Build the conda environments
RUN mamba env create --prefix /conda-envs/065992a763aefeb76352322080372458 --file /conda-envs/065992a763aefeb76352322080372458/environment.yaml && \
    mamba env create --prefix /conda-envs/5b69ffa299c6700af906ec1f423eda55 --file /conda-envs/5b69ffa299c6700af906ec1f423eda55/environment.yaml && \
    mamba env create --prefix /conda-envs/42bbaac521b8a2ad3736d2332ea6c52b --file /conda-envs/42bbaac521b8a2ad3736d2332ea6c52b/environment.yaml && \
    mamba env create --prefix /conda-envs/a5d064b61a675a69875701d952d94dbc --file /conda-envs/a5d064b61a675a69875701d952d94dbc/environment.yaml && \
    mamba clean --all -y

    # Do this LAST so code changes don't trigger a re-build of all conda envs
COPY . /app

