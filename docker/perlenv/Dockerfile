FROM nfcore/base
LABEL authors="Felipe Almeida"
LABEL description="Docker image containing perl-based bacannot tools"
ENV IMAGE=perl

# Install the conda environment
RUN conda install \
	-c bioconda -c defaults -c conda-forge -c anaconda --force-reinstall --update-deps --no-channel-priority \
	'prokka>=1.14' \
	'hmmer=3.1b2' \
	barrnap \
	mlst \
	'islandpath>=1.0.6' \
	'python>3' \
	'perl>=5.26' \
	'perl-bioperl>=1.7.8' \
	perl-digest-sha1 \
	perl-app-cpanminus \
	perl-local-lib \
	circos
ENV PERL5LIB="/opt/conda/lib/perl5/site_perl"

# Fix perl -- prokka
RUN apt-get update -y && apt-get install -y build-essential libexpat1-dev
RUN cpanm Test::Needs --force --reinstall
RUN cpanm Test::RequiresInternet Test::NoWarnings  --force --reinstall
RUN conda install -c conda-forge -y 'perl-xml-parser>2.44'
RUN cpanm XML::Twig Bio::Perl --force --reinstall || true

# set CONDA_PREFIX
ENV CONDA_PREFIX=/opt/conda

WORKDIR /work

# fix permissions for singularity
RUN chmod -R 777 /work /opt/conda/db