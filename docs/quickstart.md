# Quickstart

For a rapid and simple quickstart that enables to understand most of the available features we will use as input the _Escherichia coli_ reference genome.

## Required inputs

To run the pipeline, we basically need a samplesheet describing the genomes to be samples to be analysed (`--input`) and the path to the directory containing the databases used by bacannot (`--bacannot_db`).

## Downloading/Generating the inputs

### Input genome and samplesheet

First we need to download the genome:

```bash
# Download the ecoli ref genome
wget -O ecoli_ref.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz
gzip -d ecoli_ref.fna.gz
```

After downloading it, we must create a samplesheet for the input data as described in the [samplesheet manual page](samplesheet.md#). A proper formated file for this data would look like that:

```yaml
samplesheet: # this header is required
  - id: ecoli
    assembly: ecoli_ref.fna
    resfinder: Escherichia coli
```

!!! tip

    Download this file and save it as `bacannot_samplesheet.yaml` to help on later reference to it

### Bacannot databases

Bacannot databases are not inside the docker images anymore to avoid huge images and problems with connections and limit rates with dockerhub.

#### Pre-formatted

Users can directly download pre-formatted databases from Zenodo: https://doi.org/10.5281/zenodo.7615811

Useful for standardization and also overcoming known issues that may arise when formatting databases with `singularity` profile.

A module to download the latest pre-formatted database has also been made available:

```bash
# Download pipeline pre-built databases
nextflow run fmalmeida/bacannot --get_zenodo_db --output ./ -profile <docker/singularity>
```

#### I want to generate a new formatted database

```{bash .annotate hl_lines="5"}
# Download pipeline databases
nextflow run fmalmeida/bacannot \
    --get_dbs \
    --output bacannot_dbs \
    -profile docker
```

!!! important "About profiles"
    
    Users **must** select one of the available profiles: docker or singularity. Conda may come in future. Please read more about how to [proper select NF profiles](profiles.md#)

## Run the pipeline

In this step we will get a major overview of the main pipeline's steps. To run it, we will use the databases (`bacannot_dbs`) downloaded in the previous step.

```bash
# Run the pipeline using the Escherichia coli resfinder database
nextflow run fmalmeida/bacannot \
    --input bacannot_samplesheet.yaml \
    --output _ANNOTATION \
    --bacannot_db ./bacannot_dbs \
    --max_cpus 10 \
    -profile docker
```

!!! note "About resfinder"

    The resfinder species could also be selected via the command line with `--resfinder_species`. Please, read more about it at [manual](manual.md#) and [samplesheet](samplesheet.md#) reference pages.

### Outputs

A glimpse over the main outputs produced by bacannot is given at [outputs](outputs.md#) section.

### Testing more workflows

Moreover, we have also made available a few example datasets in the pipeline so users can test all capabilities at once, from assembling raw reads to annotating genomes. To test it users must run:

```bash
# Run the pipeline using the provided (bigger) test dataset
nextflow run fmalmeida/bacannot -profile docker,test --bacannot_db ./bacannot_dbs --max_cpus 10

# Or run the quick test
nextflow run fmalmeida/bacannot -profile docker,quicktest --bacannot_db ./bacannot_dbs ---max_cpus 10
```

!!! info ""

    Unfortunately, due to file sizes, we could not provide fast5 files for users to check on the methylation step.

### Annotation with bakta

User can also perform the core generic annotation with bakta instead of prokka. Please read [the manual](manual.md#bakta-annotation).
