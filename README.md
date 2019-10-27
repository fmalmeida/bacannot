# Bacterial Annotation (bacannot) Pipeline

This is an easy to use pipeline that uses state-of-the-art software for prokaryotic genome annotation and has only two dependencies: [Docker](https://www.docker.com/) and [Nextflow](https://github.com/nextflow-io/nextflow). Bacannot pipeline is a nextflow docker-based wrapper around a several tools that enables a better understanding of prokaryotic genomes. It uses [Prokka](https://github.com/tseemann/prokka) for generec annotation, [barrnap](https://github.com/tseemann/barrnap) for rRNA prediction. [mlst](https://github.com/tseemann/mlst) for classification within multilocus sequence types, [KofamScan](https://github.com/takaram/kofam_scan) for KO annotation, [Nanopolish](https://github.com/jts/nanopolish) for methylation annotation, [DIAMOND](https://github.com/bbuchfink/diamond) for sequence similarity searches, [JBrowse](http://jbrowse.org/) for genome browser production, [bedtools](https://bedtools.readthedocs.io/en/latest/) for gene merge, [AMRFinderPlus](https://github.com/ncbi/amr/wiki) for antimicrobial genes annotation, [Phigaro](https://github.com/bobeobibo/phigaro) and [VirSorter](https://github.com/simroux/VirSorter) for prophage sequences prediction.

<p align="center">
<img src="annotation_en.png">
<p align="center">Pipeline Fluxogram<p align="center">
</p>

## Table of contents

* [Requirements](https://github.com/fmalmeida/ngs-preprocess#requirements)
* [Quickstart](https://github.com/fmalmeida/ngs-preprocess#quickstart)
* [Documentation](https://github.com/fmalmeida/ngs-preprocess#documentation)
  * [Full usage](https://github.com/fmalmeida/ngs-preprocess#usage)
  * [Usage Examples](https://github.com/fmalmeida/ngs-preprocess#usage-examples)
  * [Configuration File](https://github.com/fmalmeida/ngs-preprocess#using-the-configuration-file)

## Requirements

* Unix-like operating system (Linux, macOS, etc)
* Java 8
* Docker
  * `fmalmeida/compgen:BACANNOT`
  * `fmalmeida/compgen:KOFAMSCAN`
  * `fmalmeida/compgen:RENV`
  * `fmalmeida/compgen:VIRSORTER`

This images have been kept separate to not create massive Docker image and to avoid dependencies conflits.

## Quickstart

1. If you don't have it already install Docker in your computer. Read more [here](https://docs.docker.com/).
    * You can give this [in-house script](https://github.com/fmalmeida/bioinfo/blob/master/dockerfiles/docker_install.sh) a try.
    * After installed, you need to download the required Docker images

          nextflow fmalmeida/bacannot --download_docker_images

2. Install Nextflow (version 0.24.x or higher):

       curl -s https://get.nextflow.io | bash

3. Give it a try:

       nextflow fmalmeida/bacannot --help

## Documentation

### Usage

    Usage:
    nextflow run fmalmeida/bacannot [--help] [ -c nextflow.config ] [OPTIONS] [-with-report] [-with-trace] [-with-timeline]

    Comments:
    This pipeline contains a massive amount of configuration variables and its usage as CLI parameters would
    cause the command to be huge.

    Therefore, it is extremely recommended to use the nextflow.config configuration file in order to make
    parameterization easier and more readable.

    Creating a configuration file:
    nextflow run fmalmeida/bacannot [--get_config]

    Show command line examples:
    nextflow run fmalmeida/bacannot --examples

    Execution Reports:
    nextflow run fmalmeida/bacannot [ -c nextflow.config ] -with-report
    nextflow run fmalmeida/bacannot [ -c nextflow.config ] -with-trace
    nextflow run fmalmeida/bacannot [ -c nextflow.config ] -with-timeline

    OBS: These reports can also be enabled through the configuration file.

    OPTIONS:

             General Parameters - Mandatory

     --outDir <string>                      Output directory name
     --threads <int>                        Number of threads to use
     --genome <string>                      Query Genome file
     --bedtools_merge_distance              Minimum number of overlapping bases for gene merge
                                            using bedtools merge.

             Prokka complementary parameters

     --prokka_center <string>               Your institude acronym to be used by prokka when
                                            renaming contigs.
     --prokka_kingdom <string>              Prokka annotation mode. Possibilities (default 'Bacteria'):
                                            Archaea|Bacteria|Mitochondria|Viruses
     --prokka_genetic_code <int>            Genetic Translation code. Must be set if kingdom is not
                                            default (in blank).
     --prokka_use_rnammer                   Tells prokka wheter to use rnammer instead of barrnap.
     --prokka_genus <string>                Set only if you want to search only a specific genus database
 
             Diamond (blastx) search parameters

     --diamond_virulence_identity           Min. identity % for virulence annotation
     --diamond_virulence_queryCoverage      Min. query coverage for virulence annotation
     --diamond_MGEs_identity                Min. identity % for ICEs and prophage annotation
     --diamond_MGEs_queryCoverage           Min. query coverage for ICEs and prophage annotation
     --diamond_minimum_alignment_length     Min. alignment length for diamond annotation

             Configure Optional processes

     --virulence_search                     Tells wheter you want or not to execute virulence annotation
     --vfdb_search                          Tells wheter you want or not to used VFDB database for virulence
                                            annotation. It is useless if virulence_search is not true
     --victors_search                       Tells wheter you want or not to used victors database for virulence
                                            annotation. It is useless if virulence_search is not true
     --resistance_search                    Tells wheter you want or not to execute resistance annotation
     --ice_search                           Tells wheter you want or not to execute ICE annotation
     --prophage_search                      Tells wheter you want or not to execute prophage annotation
     --execute_kofamscan                    Tells wheter you want or not to execute KO annotation with kofamscan

             Configure optional Methylation annotation with nanopolish
             If left blank, it will not be executed. And, with both parameters are set
             it will automatically execute nanopolish to call methylation

     --nanopolish_fast5_dir <string>         Path to directory containing FAST5 files
     --nanopolish_fastq_reads <string>       Path to fastq files (file related to FAST5 files above)

### Usage examples:

> Illumina paired end reads. Since it will always be a pattern match, example "illumina/SRR9847694_{1,2}.fastq.gz", it MUST ALWAYS be double quoted as the example below.

    ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir outputs/illumina_paired --run_shortreads_pipeline --shortreads \
    "illumina/SRR9847694_{1,2}.fastq.gz" --reads_size 2 --lighter_execute --lighter_genomeSize 4600000 --clip_r1 5 --three_prime_clip_r1 5 \
    --clip_r2 5 --three_prime_clip_r2 5 --quality_trim 30 --flash_execute

## Using the configuration file

All the parameters showed above can be, and are advised to be, set through the configuration file. When a configuration file is set the pipeline is run by simply executing `nextflow run fmalmeida/bacannot -c ./configuration-file`

Your configuration file is what will tell to the pipeline the type of data you have, and which processes to execute. Therefore, it needs to be correctly set up.

Create a configuration file in your working directory:

      nextflow run fmalmeida/ngs-preprocess --get_config

# Citation

Cite this tool as:

      Felipe Marques de Almeida. (2019, September 19). fmalmeida/ngs-preprocess: A pipeline for preprocessing NGS data from multiple sequencing platforms (Version V1.0). Zenodo. http://doi.org/10.5281/zenodo.3451406

[Prokka](https://github.com/tseemann/prokka) for generec annotation, [barrnap](https://github.com/tseemann/barrnap) for rRNA prediction. [mlst](https://github.com/tseemann/mlst) for classification within multilocus sequence types, [KofamScan](https://github.com/takaram/kofam_scan) for KO annotation, [Nanopolish](https://github.com/jts/nanopolish) for methylation annotation, [DIAMOND](https://github.com/bbuchfink/diamond) for sequence similarity searches, [JBrowse](http://jbrowse.org/) for genome browser production, [bedtools](https://bedtools.readthedocs.io/en/latest/) for gene merge, [AMRFinderPlus](https://github.com/ncbi/amr/wiki) for antimicrobial genes annotation, [Phigaro](https://github.com/bobeobibo/phigaro) and [VirSorter](https://github.com/simroux/VirSorter) for prophage sequences prediction.
