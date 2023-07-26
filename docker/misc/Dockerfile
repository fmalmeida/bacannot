FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing any-based bacannot tools"

# install mamba
RUN conda update -n root conda --yes && \
      conda update -n base conda --yes && \
      conda install -c conda-forge -y 'mamba>=1.4'

# Install the conda environment
RUN mamba install -y \
      -c bioconda -c defaults -c conda-forge -c anaconda -c falmeida \
      --no-channel-priority \
      'conda-forge::python>=3.7' \
      'blast>=2.12' \
      'diamond>=2.0.15' \
      'bedtools>=2.30' \
      'kma' \
      'nanopolish' \
      'biopython==1.78' \
      seqkit \
      bioawk \
      'easy_circos>=0.3' \
      'falmeida-py>=1.2.2' \
      'conda-forge::openssl>=1.1.1' \
      emboss \
      libtiff \
      jq && \
      mamba clean -afy

# Install samtools
RUN apt-get update -y && apt-get install -y samtools

# Install gff-toolbox
RUN git clone https://github.com/fmalmeida/gff-toolbox.git && \
      cd gff-toolbox && \
      python setup.py install && \
      gff-toolbox -h

# Create env for digIS
RUN mamba create -y \
      -c bioconda -c defaults -c conda-forge -c anaconda -c falmeida \
      -n digIS \
      --no-channel-priority \
      'hmmer==3.1b2' 'biopython==1.77' nomkl && \
	mamba clean -afy

# Install pip packages
# RUN pip install docopt pandas tabulate numpy bcbio-gff cgecore gitpython setuptools python-dateutil 'biopython==1.78'

# Install KEGGDecoder
RUN mamba create \
      -n KEGGDecoder \
      python=3.6 && \
      conda run -n KEGGDecoder python3 -m pip install KEGGDecoder && \
	mamba clean -afy

# set CONDA_PREFIX
ENV CONDA_PREFIX=/opt/conda

# get a copy of argminer database
WORKDIR /work
COPY argminer_bkp/argminer.fasta /work/argminer.fasta
COPY victors_bkp/victors_06-2022.fasta /work/victors.fasta

# get a copy of resfinder
RUN mamba create -y -n resfinder \
      -c bioconda -c defaults -c conda-forge -c anaconda -c falmeida \
      'resfinder>=4.1' docopt pandas && \
	mamba clean -afy

# get a copy of digis
RUN git clone -b master https://github.com/janka2012/digIS.git
COPY custom_fix_grange_digis.py /work/digIS/src/common/grange.py
ENV PATH=/work/digIS:$PATH

# Create env for antismash
RUN conda create -y -n antismash -c bioconda -c conda-forge \
      'antismash>=6' 'anaconda::jinja2' 'anaconda::markupsafe' nomkl && \
      rm -rf /opt/conda/envs/antismash/lib/*/site-packages/antismash/databases && \
	mamba clean -afy

# fix bioperl
RUN mamba create -n perl -y -c bioconda -c conda-forge perl-bioperl perl-app-cpanminus perl-yaml
RUN conda run -n perl PERL5LIB= PERL_LOCAL_LIB_ROOT= cpanm Bio::Root::RootI

# fix python
RUN python3 -m pip install cryptography==38.0.4

# install get zenodo
RUN pip3 install zenodo_get

# fix permissions
RUN chmod 777 -R /work
RUN chmod 777 -R /opt/conda/envs/antismash/lib/*/site-packages/antismash
RUN chmod 777 -R /opt/conda/envs/resfinder