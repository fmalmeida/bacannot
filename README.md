# Bacterial Annotation (bacannot) Pipeline

[![DOI](https://zenodo.org/badge/217119558.svg)](https://zenodo.org/badge/latestdoi/217119558) ![](https://img.shields.io/github/v/release/fmalmeida/bacannot) [![Build Status](https://travis-ci.com/fmalmeida/bacannot.svg?branch=master)](https://travis-ci.com/fmalmeida/bacannot) ![](https://img.shields.io/docker/cloud/build/fmalmeida/bacannot) [![Documentation Status](https://readthedocs.org/projects/bacannot/badge/?version=latest)](https://bacannot.readthedocs.io/en/latest/?badge=latest) ![](https://img.shields.io/badge/Nextflow-v20.01-yellowgreen)

Bacannot is an easy to use nextflow docker-based pipeline that adopts state-of-the-art software for prokaryotic genome annotation. It is a wrapper around a several tools that enables a better understanding of prokaryotic genomes. It uses:

* [Prokka](https://github.com/tseemann/prokka) for generic annotation
* [barrnap](https://github.com/tseemann/barrnap) for rRNA prediction
* [mlst](https://github.com/tseemann/mlst) for classification within multi-locus sequence types (STs)
* [KofamScan](https://github.com/takaram/kofam_scan) for KO annotation
* [KEGGDecoder](https://github.com/bjtully/BioData/tree/master/KEGGDecoder) for drawing KO annotation
* [Nanopolish](https://github.com/jts/nanopolish) for methylation annotation
* [JBrowse](http://jbrowse.org/) for genome browser production
* [bedtools](https://bedtools.readthedocs.io/en/latest/) for annotation merging
* [AMRFinderPlus](https://github.com/ncbi/amr/wiki) and [RGI](https://github.com/arpcard/rgi) for antimicrobial genes annotation
* [Phigaro](https://github.com/bobeobibo/phigaro) for prophage sequences prediction
* [IslandPath-DIMOB](https://github.com/brinkmanlab/islandpath) for genomic islands prediction
* And the databases: [CARD](https://card.mcmaster.ca/analyze/rgi), [ARGminer](https://bench.cs.vt.edu/argminer/#/classify;gene_id=A0A0Z8UZL1), [PHASTER](https://phaster.ca/), [ICEberg](https://academic.oup.com/nar/article/47/D1/D660/5165266), [Victors](http://www.phidias.us/victors/) and [VFDB](http://www.mgc.ac.cn/VFs/main.htm)

## Table of contents

* [Requirements](https://github.com/fmalmeida/bacannot#requirements)
* [Quickstart](https://github.com/fmalmeida/bacannot#quickstart)
* [Documentation](https://github.com/fmalmeida/bacannot#documentation)
  * [Full usage](https://github.com/fmalmeida/bacannot#usage)
  * [Usage Examples](https://github.com/fmalmeida/bacannot#usage-examples)
  * [Configuration File](https://github.com/fmalmeida/bacannot#using-the-configuration-file)

## Requirements

* Unix-like operating system (Linux, macOS, etc)
* Java 8
* Docker
  * `fmalmeida/bacannot:{latest, kofamscan, jbrowse, renv}`

This images have been kept separate to not create massive Docker image and to avoid dependencies conflits.

## Quickstart

1. If you don't have it already install Docker in your computer. Read more [here](https://docs.docker.com/).
    * You can give this [in-house script](https://github.com/fmalmeida/bioinfo/blob/master/dockerfiles/docker_install.sh) a try.
    * After installed, you need to download the required Docker images

          docker pull fmalmeida/bacannot:latest
          docker pull fmalmeida/bacannot:kofamscan
          docker pull fmalmeida/bacannot:jbrowse
          docker pull fmalmeida/bacannot:renv

2. Install Nextflow (version 20.07 or higher):

       curl -s https://get.nextflow.io | bash

3. Give it a try:

       nextflow fmalmeida/bacannot --help

## Documentation

### Usage

Checkout the full usage help with `nextflow run fmalmeida/bacannot --help`

Please take a time to read the [docs](https://bacannot.readthedocs.io/en/latest/?badge=latest).

### Usage examples:

> Simple annotation example through cli

    ./nextflow run main.nf --outdir TESTE --threads 3 --genome assembly.fasta --bedtools_merge_distance -20 --not_run_kofamscan

> Running with a configuration file

    ./nextflow run fmalmeida/bacannot -c nextflow.config

## Using the configuration file

All the parameters showed above can be, and are advised to be, set through the configuration file. When a configuration file is set the pipeline is run by simply executing `nextflow run fmalmeida/bacannot -c ./configuration-file`

Your configuration file is what will tell to the pipeline the type of data you have, and which processes to execute. Therefore, it needs to be correctly set up.

Create a configuration file in your working directory:

      nextflow run fmalmeida/bacannot --get_config

# Citation

    Felipe Marques de Almeida. (2020, January 25). fmalmeida/bacannot: fmalmeida/bacannot: A pipeline for an easy but comprehensive annotation of prokaryotic genomes (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.3627670

If using this tool, remember to cite the following software:

* [Prokka](https://github.com/tseemann/prokka) for generic annotation
* [barrnap](https://github.com/tseemann/barrnap) for rRNA prediction
* [mlst](https://github.com/tseemann/mlst) for classification within multi-locus sequence types (STs)
* [KofamScan](https://github.com/takaram/kofam_scan) for KO annotation
* [KEGGDecoder](https://github.com/bjtully/BioData/tree/master/KEGGDecoder) for drawing KO annotation
* [Nanopolish](https://github.com/jts/nanopolish) for methylation annotation
* [JBrowse](http://jbrowse.org/) for genome browser production
* [bedtools](https://bedtools.readthedocs.io/en/latest/) for annotation merging
* [AMRFinderPlus](https://github.com/ncbi/amr/wiki) and [RGI](https://github.com/arpcard/rgi) for antimicrobial genes annotation
* [Phigaro](https://github.com/bobeobibo/phigaro) for prophage sequences prediction
* [IslandPath-DIMOB](https://github.com/brinkmanlab/islandpath) for genomic islands prediction
* And the databases: [CARD](https://card.mcmaster.ca/analyze/rgi), [ARGminer](https://bench.cs.vt.edu/argminer/#/classify;gene_id=A0A0Z8UZL1), [PHASTER](https://phaster.ca/), [ICEberg](https://academic.oup.com/nar/article/47/D1/D660/5165266), [Victors](http://www.phidias.us/victors/) and [VFDB](http://www.mgc.ac.cn/VFs/main.htm)
