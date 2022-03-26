FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing any-based bacannot tools"

# Install the conda environment
RUN conda install -c bioconda -c defaults -c conda-forge -c anaconda -c falmeida --force-reinstall --update-deps --no-channel-priority \
            curl \
            git \
            'python=3.7' \
            'blast=2.12.0' \
            'diamond=2.0.13' \
            'bedtools=2.30' \
            'samtools=1.14' \
            'kma' \
            'kofamscan' \
            'ncbi-amrfinderplus' \
            'nanopolish' \
            'biopython==1.78' \
            gff-toolbox \
            seqkit && \
      conda clean -afy

# Create env for digIS
RUN conda create -y -n digIS -c bioconda -c defaults -c conda-forge -c anaconda -c falmeida --no-channel-priority 'hmmer==3.1b2' 'biopython==1.77' nomkl && \
		conda clean -afy

# Create env for antismash
RUN conda create -y -n antismash -c bioconda -c defaults -c conda-forge -c anaconda -c falmeida --no-channel-priority 'antismash>=6' 'anaconda::jinja2' 'anaconda::markupsafe' emboss nomkl && \
            rm -r /opt/conda/envs/antismash/lib/*/site-packages/antismash/databases && \
		conda clean -afy

# Install pip packages
RUN pip install docopt pandas KEGGDecoder tabulate numpy bcbio-gff cgecore gitpython setuptools python-dateutil 'biopython==1.78'

# set CONDA_PREFIX
ENV CONDA_PREFIX=/opt/conda

# get a copy of argminer database
WORKDIR /work
COPY argminer_bkp/argminer.fasta /work/argminer.fasta

# get a copy of resfinder
RUN git clone -b master https://git@bitbucket.org/genomicepidemiology/resfinder.git && \
      chmod a+rwx -R resfinder
ENV PATH=/work/resfinder:$PATH

# get a copy of digis
RUN git clone -b master https://github.com/janka2012/digIS.git
COPY custom_fix_grange_digis.py /work/digIS/src/common/grange.py
ENV PATH=/work/digIS:$PATH

# fix permissions
RUN chmod 777 -R /work
RUN chmod 777 -R /opt/conda/envs/antismash/lib/*/site-packages/antismash