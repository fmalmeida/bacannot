# Welcome to <u>bacannot</u> pipeline documentation

<img src="./lab_logo.png" width="300px">

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.3627669-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.3627669)
[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/fmalmeida/bacannot?include_prereleases&label=Latest%20release)](https://github.com/fmalmeida/bacannot/releases)
[![Documentation](https://img.shields.io/badge/Documentation-readthedocs-brightgreen)](https://bacannot.readthedocs.io/en/latest/?badge=latest)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![License](https://img.shields.io/badge/License-GPL%203-black)](https://github.com/fmalmeida/bacannot/blob/master/LICENSE)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40fmarquesalmeida-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/fmarquesalmeida)

## About

[Bacannot](https://github.com/fmalmeida/bacannot) is a pipeline designed to provide an easy-to-use framework for performing a comprehensive annotation on prokaryotic genomes. It is developed with [Nextflow](https://www.nextflow.io/docs/latest/index.html) and [Docker](https://www.docker.com/). It can annotate resistance genes, virulence factors, genomic islands, prophages, methylation and more.

## Workflow

The pipeline's main steps are:

| Analysis steps | Used software or databases |
| :------------- | :------------------------- |
| Genome assembly (if raw reads are given) | [Flye](https://github.com/fenderglass/Flye) and [Unicycler](https://github.com/rrwick/Unicycler) |
| Identification of closest 10 NCBI Refseq genomes | [RefSeq Masher](https://github.com/phac-nml/refseq_masher) |
| Generic annotation and gene prediction | [Prokka](https://github.com/tseemann/prokka) or [Bakta](https://github.com/oschwengers/bakta) |
| rRNA prediction | [barrnap](https://github.com/tseemann/barrnap) |
| Classification within multi-locus sequence types (STs) | [mlst](https://github.com/tseemann/mlst) |
| KEGG KO annotation and visualization | [KofamScan](https://github.com/takaram/kofam_scan) and [KEGGDecoder](https://github.com/bjtully/BioData/tree/master/KEGGDecoder) |
| Annotation of secondary metabolites | [antiSMASH](https://docs.antismash.secondarymetabolites.org/) |
| Methylation annotation | [Nanopolish](https://github.com/jts/nanopolish) |
| Annotation of antimicrobial (AMR) genes | [AMRFinderPlus](https://github.com/ncbi/amr/wiki), [ARGminer](https://bench.cs.vt.edu/argminer), [Resfinder](https://cge.cbs.dtu.dk/services/ResFinder/) and [RGI](https://github.com/arpcard/rgi) |
| Annotation of virulence genes | [Victors](http://www.phidias.us/victors/) and [VFDB](http://www.mgc.ac.cn/VFs/main.htm) |
| Prophage sequences and genes annotation | [PHASTER](http://phast.wishartlab.com/), [Phigaro](https://github.com/bobeobibo/phigaro) and [PhySpy](https://github.com/linsalrob/PhiSpy) |
| Annotation of integrative and conjugative elements | [ICEberg](http://db-mml.sjtu.edu.cn/ICEberg/) |
| Focused detection of insertion sequences | [digIS](https://github.com/janka2012/digIS) |
| _In silico_ detection of plasmids | [Plasmidfinder](https://cge.cbs.dtu.dk/services/PlasmidFinder/) and [Platon](https://github.com/oschwengers/platon) |
| Prediction and visualization of genomic islands | [IslandPath-DIMOB](https://github.com/brinkmanlab/islandpath) and [gff-toolbox](https://github.com/fmalmeida/gff-toolbox) |
| Custom annotation from formatted FASTA or NCBI protein IDs | [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs) |
| Merge of annotation results | [bedtools](https://bedtools.readthedocs.io/en/latest/) |
| Genome Browser renderization | [JBrowse](http://jbrowse.org/) |
| Circos plot generation | [easy_circos](https://easy_circos.readthedocs.io/en/latest/index.html) |
| Renderization of automatic reports and shiny app for results interrogation | [R Markdown](https://rmarkdown.rstudio.com/), [Shiny](https://shiny.rstudio.com/) and [SequenceServer](https://sequenceserver.com/) |

!!! note "Quickstart"

    A [quickstart](quickstart.md#) is available so you can quickly get the gist of the pipeline's capabilities.

!!! info "About prokka annotation"

    In order to increase the accuracy of prokka annotation, this pipeline includes an additional HMM database to prokka's defaults. It can be either TIGRFAM (smaller but curated) or PGAP (bigger comprehensive NCBI database that contains TIGRFAM).

## Usage

The pipeline's common usage is very simple as shown below:

```bash
# usual command-line
nextflow run fmalmeida/bacannot \
    --bacannot_db "./bacannot_databases" \
    --input "bacannot_samplesheet.yml"
```

!!! quote

    Some parameters are required, some are not. Please read the pipeline's manual reference to understand each parameter.

## Support contact

Whenever a doubt arise feel free to contact me at almeidafmarques@gmail.com
