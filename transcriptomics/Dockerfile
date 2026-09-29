FROM condaforge/miniforge3:latest

LABEL stage="6701-rnaseq"

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml -n 6701-rnaseq && \
    conda clean -afy

ENV PATH="/opt/conda/envs/6701-rnaseq/bin:$PATH"

# avoid sra-tools interactive config prompt in container
RUN mkdir -p /root/.ncbi && \
    printf '/LIBS/GUID = "%s"\n/repository/user/main/public/root = "/data/sra"\n' "$(cat /proc/sys/kernel/random/uuid)" > /root/.ncbi/user-settings.mkfg

WORKDIR /app