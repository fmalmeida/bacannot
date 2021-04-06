# bacannot pipeline changelog

The tracking for changes started in v2.1

## v2.2.1

Added a parameter `--singularity` so that users can easily use Singularity instead of Docker, if wanted.

## v2.2

### Automatic reports

Some problems with the automatic reports, such as not understanding missing or empty files, and conditionals between processes, have been addressed and fixed.

### VFDB database

VFDB database has been set again to use only the core (A) virulence gene dataset. The use of complete (B) dataset were returning lots of spurious hits.

### Resistance annotation

The Resfinder software and the ARGminer database have been added to the pipeline in order to augment the diversity of sources for resistance genes annotation.

### Prophage annotation

The PhiSpy software has been added to the pipeline to increase the diversity of sources for prophage sequences annotation. This software is complementary to Phigaro.

### KO annotation

KeggDecoder has been added to provide a quick and simple visualization of the KEGG Orthology (KO) annotation results.

### Plasmid annotation

The Platon software has been added to aid in the prediction of plasmid sequences, together with plasmidfinder.

### Custom annotation

The possibility to perform additional annotations using user's custom nucleotide gene databases has been added with the `--custom_db` parameter. The proper configuration of these databases are documented at: https://bacannot.readthedocs.io/en/latest/custom-db.html

### Multiple genome analysis

The possibility to perform the annotation of multiple genomes at once has been added with the `--in_yaml` parameter which takes as input a samplesheet file in YAML format as described in the documentation at: https://bacannot.readthedocs.io/en/latest/samplesheet.html

### Annotation from raw reads

The possibility to annotate genomes from raw reads have been added with the parameters `--sreads_single`, `--sreads_paired` and `--lreads` which takes as input raw sequencing reads from Illumina, Pacbio and ONT in FASTq format and uses it to assemble the genome with Unicycler or Flye (depending on the data availability) prior to the annotation step.

### Bacannot shiny server

A simple shiny server has been created and implemented with the `run_server.sh` bash script that loads the shiny app from a docker image that allows the user to quickly interrogate the annotation results via the JBrowse genome browser, the annotation reports and with a BLAST tool implemented in the server that enables users to quickly detect the presence of additional genes/sequences. Take a better look at: https://bacannot.readthedocs.io/en/latest/outputs.html

## v2.1

This versions have a few additions to the pipeline workflow, they are highlighted and explained below:

### nf-core schema

We have added a nextflow parameter schema in json that is compliant with nf-core. This enables that users trigger the graphical interface for configuration and execution of the pipeline via [nf-core launch](https://nf-co.re/launch) utility, also it is possible to track the pipeline execution with [nextflow tower](https://tower.nf/).

```bash
# It is triggered as
nf-core launch fmalmeida/bacannot
```

### plasmidfinder

_In silico_ plasmid detection has been added with Plasmidfinder. In order to not execute this parameter, users will need to use `--not_run_plasmid_search`. Otherwise, its thresholds can be configured with `--plasmids_minid` and `--plasmids_mincov`.

### VFDB database

VFDB database has been changed from the core (A) dataset to the complete (B) dataset. This dataset is first clustered with cd-hit using a 90% identity threshold and this new database file (after cd-hit) is used for the virulence gene annotation.
