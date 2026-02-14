FROM condaforge/miniforge3:latest
WORKDIR /tmp
SHELL ["/bin/bash", "-l", "-c"]

# Install minimal dependencies (Snakemake)
RUN mamba install -c bioconda -c conda-forge snakemake -y

# Copy workflow files
COPY ./workflow /tmp/workflow
COPY ./config /tmp/config

# Note: The 'src' directory (containing external scripts like remove_rec...) should be mounted by the user
# to /tmp/program if needed.

CMD [ "/bin/bash" ]