FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing python-based bacannot tools"

# Install the conda environment
# RUN conda install -c conda-forge -y 'mamba>=0.21'
RUN conda install -y \
      -c bioconda -c defaults -c conda-forge -c anaconda \
      'openssl=1.1.1' \
      'platon>=1.6' \
      phispy \
      plasmidfinder \
      'python>=3.7' \
      refseq_masher
RUN conda install -c conda-forge -y 'mamba>=0.21'
RUN mamba create -n rgi -c bioconda -c defaults -c conda-forge -c anaconda 'rgi>=5.2.1'
RUN mamba create -n phigaro -c bioconda -c defaults -c conda-forge -c anaconda phigaro

# set CONDA_PREFIX
ENV CONDA_PREFIX=/opt/conda
