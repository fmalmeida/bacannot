FROM ubuntu:22.04

LABEL MAINTAINER="Felipe Marques de Almeida <marques.felipe@aluno.unb.br>"

# Workdir
WORKDIR /work

# Set main R dependencies
## Add R repo key
RUN apt-get update && \
			DEBIAN_FRONTEND=noninteractive apt-get install -y software-properties-common gnupg build-essential bedtools && \
			apt-get update

## Install R
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y tzdata && \
			DEBIAN_FRONTEND=noninteractive apt-get update && apt-get install -y r-base r-base-core

## Install R-packages
RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
			DEBIAN_FRONTEND=noninteractive apt-get --option Acquire::Retries=100 --option Acquire::http::Timeout="300" install -yq r-cran-optparse r-cran-rmysql r-bioc-gviz r-cran-xml pandoc
RUN DEBIAN_FRONTEND=noninteractive apt-get --option Acquire::Retries=100 --option Acquire::http::Timeout="300" install -yq r-cran-plyr libssl-dev libcurl4-openssl-dev libxml2-dev libudunits2-dev libcairo2-dev libpq-dev r-cran-rmarkdown
ADD scripts/rscripts/installPack.R /work/rscripts/installPack.R
RUN DEBIAN_FRONTEND=noninteractive Rscript /work/rscripts/installPack.R && \
			apt-get clean

# METADATA and REPORTS
# Add local files
ADD reports /work/reports

# Download Victors metadata
RUN apt-get install -y wget && \
			cd /work/reports && \
			( wget "http://www.phidias.us/victors/export_text.php?c_mc_pathogen_id=&c_phi_function=Virulence%20factor&c_mc_victor_name=&c_gene_locus_tag=&db_type=gene_gi&db_id=&c_max_tmhmm_PredHel=1&c_max_tmhmm_PredHel_check=&c_min_spaan_score_check=&c_min_spaan_score=0.51&keywords=&c_human_alignment=&c_mouse_alignment=&c_localization[]=Any&cog_cat_id[]=" -O /work/reports/new_victors_metadata.tsv && cat /work/reports/new_victors_metadata.tsv > /work/reports/victors_metadata.tsv ) || ( echo "using backup" )

# Download CARD metadata
RUN cd /work/reports && \
			wget -O card-data.tar.bz2 https://card.mcmaster.ca/download/0/broadstreet-v3.1.0.tar.bz2 && \
			tar jxvf card-data.tar.bz2 && \
			rm card-data.tar.bz2

# Get scripts
ADD scripts/rscripts /work/rscripts
ADD scripts/bscripts /work/bscripts

# Make Rmds executable
RUN chmod a+rx /work/reports/*.Rmd

# Make rscripts executable
RUN chmod a+rx /work/rscripts/*

# Make bscripts executable
RUN chmod a+rx /work/bscripts/*

# Set PATH
ENV PATH="$PATH:/work/rscripts"

# Re-set workidir
WORKDIR /work
RUN chmod -R a+rw /work
