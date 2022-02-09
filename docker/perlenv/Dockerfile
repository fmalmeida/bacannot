FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing perl-based bacannot tools"
ENV IMAGE=perl
# Install the conda environment
COPY ${IMAGE}.yml /
RUN conda env create --quiet -f /${IMAGE}.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/bacannot-${IMAGE}/bin:$PATH

# Fix perl -- prokka
RUN conda update prokka -n bacannot-${IMAGE}
RUN apt-get update -y && apt-get install -y build-essential
RUN cpanm Bio::SearchIO::hmmer3 --force


# Dump the details of the installed packages to a file for posterity
RUN conda env export --name bacannot-${IMAGE} > bacannot-${IMAGE}.yml

# set CONDA_PREFIX
ENV CONDA_PREFIX=/opt/conda
