# Welcome to <u>phylogram</u> pipeline documentation

<img src="./lab_logo.png" width="300px">

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40fmarquesalmeida-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/fmarquesalmeida)

## About

[Phylogram](https://github.com/fmalmeida/phylogram) is a pipeline that automatically reconstructs phylogeny relationships and trees based on domain alignments. It aligns domains to sequences using `hmmsearch` and creates trees and annotation files ready to be loaded in [Itol](https://itol.embl.de/).

## Workflow

The steps taken by the pipeline to reconstruct the phylogeny between sequences is outlined below:

1. Subset a references from [RefPlants](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001124)
2. Aligns input sequences against desired Pfam domain with [hmmsearch](https://github.com/EddyRivasLab/hmmer)
3. Alignment is "preprocessed" with [ClipKIT](https://github.com/JLSteenwyk/ClipKIT)
4. Alignment is trimmed with [Trimal](https://github.com/inab/trimal)
5. Alignment is weighted with [Tcoffee TCS](https://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#transitive-consistency-score-tcs)
6. Phylogeny tree is built with [FastTree](http://www.microbesonline.org/fasttree/), [IQTree](http://www.iqtree.org/) or [Raxml-ng](https://github.com/amkozlov/raxml-ng)
    + Best model is selected with [ModelFinder](http://dx.doi.org/10.1038/nmeth.4285) or [jmodeltest](https://github.com/ddarriba/jmodeltest2)
7. [TreeShrink](https://github.com/uym2/TreeShrink) and [gotree](https://github.com/evolbioinfo/gotree) are used to prune the final tree, removing long branches
8. Itol Annotation files are generated with custom scripts

!!! note "Quickstart"

    A [quickstart](quickstart.md#) is available so you can quickly get the gist of the pipeline's capabilities.

## Usage

The pipeline's common usage is very simple as shown below:

```bash
# usual command-line
nextflow run fmalmeida/phylogram \
    --prefix "testing" \
    --fasta "input.faa" \
    --hmm "PF00931" \
    --refplants_species "Arabidopsis"
```

> Some parameters are required, some are not. Please read the pipeline's manual reference to understand each parameter.

## Testing

A testing profile is available with:

```bash
nextflow run fmalmeida/phylogram -profile docker,test
```
