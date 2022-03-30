FROM ubuntu:20.04

LABEL MAINTAINER Felipe Marques de Almeida <felipemarques89@gmail.com>

WORKDIR /work
RUN chmod 777 -R /work

RUN apt-get update && apt-get upgrade -y
RUN apt-get -qq update --fix-missing
RUN export DEBIAN_FRONTEND="noninteractive" ; apt-get --no-install-recommends -y install git build-essential zlib1g-dev libxml2-dev libexpat-dev postgresql-client libpq-dev libpng-dev curl wget unzip rsync

RUN wget --no-check-certificate https://github.com/GMOD/jbrowse/releases/download/1.16.9-release/JBrowse-1.16.9.zip && \
    unzip JBrowse-1.16.9.zip && \
    mv JBrowse-1.16.9 jbrowse && \
    apt-get install -y libnet-ssleay-perl samtools && \
    cd jbrowse && ./setup.sh && \
    ./bin/cpanm --force JSON Hash::Merge PerlIO::gzip Devel::Size Heap::Simple Heap::Simple::XS List::MoreUtils Exception::Class Test::Warn Bio::Perl Bio::DB::SeqFeature::Store File::Next Bio::DB::Das::Chado && \
    perl Makefile.PL && make && make install

ENV PATH=/work/jbrowse/bin:$PATH

RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig ; mv bedGraphToBigWig /work/jbrowse/bin ; chmod a+x /work/jbrowse/bin/bedGraphToBigWig

# Install FeatureIO
RUN cd /work/jbrowse && \
    ./bin/cpanm --force Bio::FeatureIO && \
    ./bin/cpanm --force Bio::GFF3::LowLevel Bio::GFF3::LowLevel::Parser DBD::SQLite \
    Digest::Crc32 File::Copy::Recursive JSON::XS Parse::RecDescent local::lib

# for backup reasons to work with previous versions
RUN mkdir -p /work/bscripts && \
      wget -O /work/bscripts/run_jbrowse.sh https://github.com/fmalmeida/bacannot/raw/305956e6ae2b42a8a5f69bedccbacccadaf2a30e/docker/scripts/bscripts/run_jbrowse.sh && \
      chmod a+x /work/bscripts/run_jbrowse.sh

# Change Permissions
RUN chmod 777 -R /work

WORKDIR /work/data
