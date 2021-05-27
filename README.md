[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3627669.svg)](https://doi.org/10.5281/zenodo.3627669) [![Releases](https://img.shields.io/github/v/release/fmalmeida/bacannot)](https://github.com/fmalmeida/bacannot/releases) [![Documentation](https://img.shields.io/badge/Documentation-readthedocs-brightgreen)](https://bacannot.readthedocs.io/en/latest/?badge=latest) [![Dockerhub](https://img.shields.io/badge/Docker-fmalmeida/bacannot-informational)](https://hub.docker.com/r/fmalmeida/bacannot) [![Docker build](https://img.shields.io/docker/cloud/build/fmalmeida/bacannot)](https://hub.docker.com/r/fmalmeida/bacannot) [![Nextflow version](https://img.shields.io/badge/Nextflow%20>=-v20.07-important)](https://www.nextflow.io/docs/latest/getstarted.html) [![License](https://img.shields.io/badge/License-GPL%203-black)](https://github.com/fmalmeida/bacannot/blob/master/LICENSE)

<p align="center">

  <h1 align="center">bacannot pipeline</h2>

  <p align="center">
    <h3 align="center">A generic but comprehensive bacterial annotation pipeline</h3>
    <br />
    <a href="https://bacannot.readthedocs.io/en/latest/index.html"><strong>See the documentation »</strong></a>
    <br />
    <br />
    <a href="https://github.com/fmalmeida/bacannot/issues">Report Bug</a>
    ·
    <a href="https://github.com/fmalmeida/bacannot/issues">Request Feature</a>
  </p>
</p>

## About

Bacannot is an easy to use nextflow docker-based pipeline that adopts state-of-the-art software for prokaryotic genome annotation. It is a wrapper around a several tools that enables a better understanding of prokaryotic genomes.

See [the list of tools](list_of_tools.md).

### Release notes

Are you curious about changes between releases? See the [changelog](CHANGELOG.md).

* If you want to execute a release that is not the latest, you must build the docker image of that specific release since available tools and image strucure may change. [See the instructions for that](earlier_releases_instructions.md).
* Versions below 2.0 are no long supported.

### Complementary pipelines

Moreover, this pipeline has two complementary pipelines (also written in nextflow) for [NGS preprocessing](https://github.com/fmalmeida/ngs-preprocess) and [Genome assembly](https://github.com/fmalmeida/MpGAP) that can give the user a complete workflow for bacterial genomics analyses.

## Requirements

* Unix-like operating system (Linux, macOS, etc)
  + Windows users maybe can execute it using the linux subsystem for windows as shown in:
    + https://docs.microsoft.com/pt-br/windows/wsl/install-win10
    + https://www.nextflow.io/docs/latest/getstarted.html
    + https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux
* Java 8
* Docker
  + `fmalmeida/bacannot:{latest, kofamscan, jbrowse, renv}`

These images have been kept separate to not create massive Docker image and to avoid dependencies conflicts.

## Installation

1. If you don't have it already install [Docker](https://docs.docker.com/) in your computer.
    * After installed, you need to download the required Docker images

          docker pull fmalmeida/bacannot:latest
          docker pull fmalmeida/bacannot:kofamscan
          docker pull fmalmeida/bacannot:jbrowse
          docker pull fmalmeida/bacannot:renv
          docker pull fmalmeida/bacannot:server (For the shiny parser)
          docker pull fmalmeida/mpgap (Only necessary if using raw reads as input)

2. Install Nextflow (version 20.07 or higher):

       curl -s https://get.nextflow.io | bash

3. Give it a try:

       nextflow run fmalmeida/bacannot --help

> Users can get let the pipeline always updated with: `nextflow pull fmalmeida/bacannot`

### Maintaining databases up-to-date

By default, github actions have been set to build the docker image containing the databases (`fmalmeida/bacannot:latest`) in the first day of every month. Therefore, to use the most up-to-date databases users must run `docker pull fmalmeida/bacannot:latest` before running the pipeline.

Additionally, a custom script is provided to allow users to update the database image any time.

```bash
bash <(wget -O - -o /dev/null https://github.com/fmalmeida/bacannot/raw/develop/bin/update_database_image.sh)
```

> This command line will trigger a custom script that downloads the databases and build the fmalmeida/bacannot:latest docker image.

## Quickstart

For a rapid and simple quickstart we will use as input the nanopore raw reads provided in the [Canu quickstart section](https://canu.readthedocs.io/en/latest/quick-start.html#assembling-pacbio-clr-or-nanopore-data).

```bash

  # Download the data and save it as oxford.fasta
  curl -L -o oxford.fasta http://nanopore.s3.climb.ac.uk/MAP006-PCR-1_2D_pass.fasta

  # Run the pipeline using the Escherichia coli resfinder database
  nextflow run fmalmeida/bacannot \
  --prefix ecoli \
  --lreads oxford.fasta \
  --lreads_type nanopore \
  --outdir _ANNOTATION \
  --threads 4 \
  --resfinder_species "Escherichia coli"

```

### Overview of outputs

A nice overview of the output directory structure and the main tools/features produced by the pipeline is provided at https://bacannot.readthedocs.io/en/latest/outputs.html.

## Documentation

### Usage

<a href="https://bacannot.readthedocs.io/en/latest/index.html"><strong>Users are advised to read the complete documentation »</strong></a>

* Complete command line explanation of parameters:
    + `nextflow run fmalmeida/bacannot --help`
* See usage examples in the command line:
    + `nextflow run fmalmeida/bacannot --examples`

### Command line usage examples

Command line executions are exemplified [in the manual](https://bacannot.readthedocs.io/en/latest/examples.html).

#### Using the configuration file

All the parameters showed above can be, and are advised to be, set through the configuration file. When a configuration file is set the pipeline is run by simply executing `nextflow run fmalmeida/bacannot -c ./configuration-file`

Your configuration file is what will tell to the pipeline the type of data you have, and which processes to execute. Therefore, it needs to be correctly set up.

Create a configuration file in your working directory:

      nextflow run fmalmeida/bacannot --get_config

### Interactive graphical configuration and execution

#### Via NF tower launchpad (good for cloud env execution)

Nextflow has an awesome feature called [NF tower](https://tower.nf). It allows that users quickly customise and set-up the execution and configuration of cloud enviroments to execute any nextflow pipeline from nf-core, github (this one included), bitbucket, etc. By having a compliant JSON schema for pipeline configuration it means that the configuration of parameters in NF tower will be easier because the system will render an input form.

Checkout more about this feature at: https://seqera.io/blog/orgs-and-launchpad/

<p align="center">
<img src="https://j.gifs.com/GRnqm7.gif" width="500px"/>
</p>

#### Via nf-core launch (good for local execution)

Users can trigger a graphical and interactive pipeline configuration and execution by using [nf-core launch](https://nf-co.re/launch) utility. nf-core launch will start an interactive form in your web browser or command line so you can configure the pipeline step by step and start the execution of the pipeline in the end.

```bash
# Install nf-core
pip install nf-core

# Launch the pipeline
nf-core launch fmalmeida/bacannot
```

It will result in the following:

<p align="center">
<img src="https://github.com/fmalmeida/bacannot/raw/master/images/nf-core-asking.png" width="500px"/>
</p>

<p align="center">
<img src="https://github.com/fmalmeida/bacannot/raw/master/images/nf-core-gui.png" width="500px"/>
</p>

## Known issues

1. Sometimes when navigating through the shiny parser the reports and JBrowse tabs may still be pointing to old, or just different, samples that have been analysed before and not the actual sample in question. For example, you open the shiny server for the Sample 2, but the reports and JBrowse are showing results of Sample 1. This is caused by the browser's data storages and cookies.
    * To solve this problem user's can just clear the cookies and data cache from the browser.
2. The JBrowse wrapper in the shiny server is not capable of displaying the GC content and methylation plots when available. It can only display the simpler tracks. If the user wants to visualise and interrogate the GC or methylation tracks it must open the JBrowse outside from the shiny server. For that, two options are available:
    * You can navigate to the `jbrowse` directory under your sample's output folder and simply execute `http-server`. This command can be found at: https://www.npmjs.com/package/http-server
    * Or, you can download the [JBrowse Desktop app](https://jbrowse.org/docs/jbrowse_desktop.html) and, from inside the app, select the folder `jbrowse/data` that is available in your sample's output directory.

## Citation

Please cite this pipeline using our Zenodo tag or directly via the github url. Also, whenever its outputs are used/helpful, **remember** to cite the correct the used software (Check [the list of tools](list_of_tools.md)).
