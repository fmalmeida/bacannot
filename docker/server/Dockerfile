# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:latest

# Owner
LABEL MAINTAINER Felipe Marques de Almeida <marques.felipe@aluno.unb.br>

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    r-bioc-biostrings \
    r-bioc-iranges

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean

## install ncbi-blast
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz && \
    tar zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz && \
    rm ncbi-blast-2.13.0+-x64-linux.tar.gz && \
    mv ./ncbi-blast-2.13.0+/bin/* /usr/local/bin/ && \
    rm -rf ./ncbi-blast-2.13.0+/

# Set workdir
WORKDIR /app

# Install Java
RUN apt-get -y install default-jre r-cran-rjava

# Install renv & restore packages
RUN Rscript -e 'install.packages(c("bookdown", "rmarkdown", "shiny", "shinyFiles", "dplyr", "dbplyr", "RSQLite", "DT", "glue", "stringr", "DataCombine", "shinythemes", "stringr", "openxlsx"), dependencies=TRUE)'

# Install bedtools
RUN apt-get install -y bedtools

# Install curl
RUN apt-get install -y curl

# Install sequenceserver
RUN apt-get install -y ruby-full
RUN gem install sequenceserver

## Get scripts
COPY scripts/shiny_parser.Rmd /app/shiny_parser.Rmd
COPY scripts/shiny_parser_func.R /app/shiny_parser_func.R
COPY scripts/get_attributes_field.R /app/get_attributes_field.R
COPY scripts/blast_app.R /app/blast_app.R
COPY scripts/SQL_parser.R /app/SQL_parser.R
COPY scripts/server_entry.sh /app/server_entry.sh

# expose port
EXPOSE 3838
EXPOSE 4567

# Prevent SequenceServer from prompting user to join announcements list.
RUN mkdir -p ~/.sequenceserver && touch ~/.sequenceserver/asked_to_join

# run app on container start
CMD bash /app/server_entry.sh
