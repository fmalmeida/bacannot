FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing perl-based bacannot tools"
ENV IMAGE=perl
# Install the conda environment
RUN conda install -y -c conda-forge mamba
COPY ${IMAGE}.yml /
RUN mamba env create --quiet -f /${IMAGE}.yml && mamba clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/bacannot-${IMAGE}-v3.1/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name bacannot-${IMAGE}-v3.1 > bacannot-${IMAGE}-v3.1.yml

# set CONDA_PREFIX
ENV CONDA_PREFIX=/opt/conda