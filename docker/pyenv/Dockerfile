FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing perl-based bacannot tools"

# Install the conda environment
RUN conda install -c bioconda -c defaults -c conda-forge -c anaconda --force-reinstall --update-deps --no-channel-priority \
      'python>3' \
	plasmidfinder \
      platon \
      phispy \
      phigaro

# set CONDA_PREFIX
ENV CONDA_PREFIX=/opt/conda
