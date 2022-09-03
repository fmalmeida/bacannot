FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing perl-based bacannot tools"
ENV IMAGE=perl

# Install the conda environment
RUN conda install -c bioconda -c defaults -c conda-forge -c anaconda --force-reinstall --update-deps --no-channel-priority \
	'prokka>=1.14' \
	'hmmer=3.1b2' \
	barrnap \
	mlst \
	islandpath \
	'python>3' \
	perl-bioperl \
	perl-bioperl-core \
	perl-digest-sha1

# Fix perl -- prokka
RUN conda update -y prokka
RUN apt-get update -y && apt-get install -y build-essential
RUN cpanm Bio::SearchIO::hmmer3 --force

# Fix perl -- bioperl
RUN cpanm Cache::FileCache --force
RUN cpanm Bio::DB::GenPept --force
RUN cpanm Bio::Perl --force

# set CONDA_PREFIX
ENV CONDA_PREFIX=/opt/conda

WORKDIR /work

# fix permissions for singularity
RUN chmod -R 777 /work /opt/conda/db