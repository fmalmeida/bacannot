FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing any-based bacannot tools"
ENV CHANNELS='-c bioconda -c defaults -c conda-forge -c anaconda -c falmeida'

# install mamba
RUN conda update -n base -c defaults conda && conda install -c conda-forge -y 'mamba>=1.0'

# Install the conda environment
RUN mamba install -y ${CHANNELS} \
      'python>=3.7' \
      'blast>=2.12' \
      'diamond>=2.0.15' \
      'bedtools>=2.30' \
      'samtools' \
      'kma' \
      'nanopolish' \
      'biopython==1.78' \
      gff-toolbox \
      seqkit \
      bioawk \
      'easy_circos>=0.3' \
      jq && \
      mamba clean -afy

# Create env for digIS
RUN mamba create -y ${CHANNELS} \
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
RUN conda create -y -n resfinder ${CHANNELS} \
      'resfinder>=4.1' docopt pandas && \
	conda clean -afy

# get a copy of digis
RUN git clone -b master https://github.com/janka2012/digIS.git
COPY custom_fix_grange_digis.py /work/digIS/src/common/grange.py
ENV PATH=/work/digIS:$PATH

# Create env for antismash
RUN mamba create -y -n antismash ${CHANNELS} \
      'antismash>=6' 'anaconda::jinja2' 'anaconda::markupsafe' emboss nomkl && \
      rm -r /opt/conda/envs/antismash/lib/*/site-packages/antismash/databases && \
	mamba clean -afy

# fix permissions
RUN chmod 777 -R /work
RUN chmod 777 -R /opt/conda/envs/antismash/lib/*/site-packages/antismash
RUN chmod 777 -R /opt/conda/envs/resfinder