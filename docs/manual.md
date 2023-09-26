# Manual

```bash
# Get help in the command line
nextflow run fmalmeida/bacannot --help
```

!!! tip

    All these parameters are configurable through a configuration file. We encourage users to use the configuration file since it will keep your execution cleaner and more readable. See a [config example](config.md#).

## Input description

### Required

To execute the annotation pipeline users **must** provide genomic data as either raw reads or assembled genomes as input. When raw reads are used, Unicycler and Flye assemblers are used to create, respectively, shortreads-only and hybrid assemblies, or longreads-only assemblies for the annotation process. Which means, the minimum required input files are:

* An assembled genome in FASTA format, **or**;
* Raw sequencing reads.

### Optional

The pipeline accepts as input two other input files types that are used to perform additional annotation processes, they are:

* path to a directory of FAST5
    * Then used together with nanopore reads it will call DNA methylation with Nanopolish.
* path to custom databases as described in [custom-db reference page](custom-db.md#)
    * These custom databases will be used to perform additional annotation processes using BLAST. Please check the both the explanation [about the parameters](manual.md#custom-nucl-databases) and about its [configuration](custom-db.md#).

## Input/output options

| <div style="width:100px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--input`  | :material-check: | NA       | Input samplesheet describing all the samples to be analysed |
| `--enable_deduplication` | :material-close: | false | Run deduplication command on input reads before assembly |
| `--output` | :material-check: | results  |  Name of directory to store output values. A sub-directory for each genome will be created inside this main directory. |
| `--bacannot_db` | :material-check: | NA | Path for root directory containing required bacannot databases |

!!! note "About the samplesheet"
    
    Please read the [samplesheet manual page](samplesheet.md#) to better understand its format.

## Database download options

| <div style="width:120px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--get_dbs`  | :material-close: | false  | Instead of running the analysis workflow, it will try to download required databases and save them in `--output` |
| `--force_update` | :material-close: | false | Instead of only downloading missing databases, download everything again and overwrite. |

!!! tip ""
    
    The quickstart shows a common usage of these parameters.

## Prokka annotation

| <div style="width:160px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--prokka_kingdom`      | :material-close: | Bacteria | Prokka annotation mode. Possibilities: Archaea|Bacteria |
| `--prokka_genetic_code` | :material-close: | 11 | Genetic Translation code. Must be set if a different kingdom is customized. |
| `--prokka_use_rnammer`  | :material-close: | false | Tells Prokka whether to use rnammer instead of barrnap |
| `--prokka_use_pgap`     | :material-close: | false | Include comprehensive PGAP hmm database in prokka annotation instead of TIGRFAM. Although comprehensive it increases runtime |

!!! info "About prokka annotation"

    In order to increase the accuracy of prokka annotation, this pipeline includes an additional HMM database to prokka's defaults. It can be either TIGRFAM (smaller but curated) or PGAP (bigger comprehensive NCBI database that contains TIGRFAM).

## Bakta annotation

!!! info "Using Bakta"

    If desired, users can use [`bakta`](https://github.com/oschwengers/bakta) instead of `prokka` to perform the core generic annotation of their prokaryotic genomes. For that, users must simply [download and store bakta database](https://github.com/oschwengers/bakta#database-download) in their machine, and pass its path to `bacannot` with the `--bakta_db` parameter.

    We opted for having it like this because bakta database is quite big.

| <div style="width:160px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--bakta_db`      | :material-close: | NA | Path to bakta database. If given, bacannot will use bakta instead of prokka. |

## Resfinder annotation

The use of this parameter sets a default value for input samples. If a sample has a different value given inside the samplesheet, the pipeline will use, for that sample, the value found inside the [samplesheet](samplesheet.md#).

| <div style="width:160px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--resfinder_species` | :material-close: | NA | Resfinder species panel. It activates the resfinder annotation process using the given species panel. Check the available species at [their main page](https://cge.cbs.dtu.dk/services/ResFinder/) and in [their repository page](https://bitbucket.org/genomicepidemiology/resfinder/src/master/#usage). If your species is not available in Resfinder panels, you may use it with the "Other" panel (`--resfinder_species "Other"`). |

## On/Off processes

| <div style="width:180px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--skip_virulence_search` | :material-close: | false | Tells whether not to run virulence factors annotation. It skips both vfdb and victors annotation |
| `--skip_plasmid_search` | :material-close: | false | Tells whether not to run plasmid detection/typing modules |
| `--skip_resistance_search` | :material-close: | false | Tells whether not to run resistance genes annotation modules |
| `--skip_iceberg_search` | :material-close: | false | Tells whether not to run mobile genetic elements annotation with ICEberg |
| `--skip_prophage_search` | :material-close: | false | Tells whether not to run prophage annotation modules |
| `--skip_kofamscan` | :material-close: | false | Tells whether not to run KEGG orthology (KO) annotation with KofamScan |
| `--skip_antismash` | :material-close: | false | Tells whether or not to run antiSMASH (secondary metabolite) annotation. AntiSMASH is executed using only its core annotation modules in order to keep it fast. |

## Custom databases

Users can give fasta files (nucl or prot) properly formatted or a text file containing a list of NCBI protein IDs (one per line). Please check the [custom db manual](custom-db.md#) for more information. Sequences are searched against the genome, with `blastn` for nucl sequences and `tblastn` for prot sequences.

| <div style="width:180px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--custom_db` | :material-close: | NA | Custom gene nucleotide/protein databases to be used for additional annotations. N files are accepted separated by commas. E.g. `--custom_db db1.fasta,db2.fasta,db3.fasta`. |
| `--ncbi_proteins` | :material-close: | NA | Path to file with NCBI protein IDs. The pipeline will download, format and use them for additional annotation. |

## Annotation thresholds

| <div style="width:200px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--blast_virulence_minid` | :material-close: | 90 | Identity (%) threshold to be used when annotating virulence factors from VFDB and Victors |
| `--blast_virulence_mincov` | :material-close: | 90 | Coverage (%) threshold to be used when annotating virulence factors from VFDB and Victors |
| `--blast_resistance_minid` | :material-close: | 90 | Identity (%) threshold to be used when annotating AMR genes with CARD-RGI, Resfinder, ARGminer and AMRFinderPlus. |
| `--blast_resistance_mincov` | :material-close: | 90 | Coverage (%) threshold to be used when annotating AMR genes with Resfinder, ARGminer and AMRFinderPlus. CARD-RGI is not affected. |
| `--plasmids_minid` | :material-close: | 90 | Identity (%) threshold to be used when detecting plasmids with Plasmidfinder |
| `--plasmids_mincov` | :material-close: | 60 | Coverage (%) threshold to be used when detecting plasmids with Plasmidfinder |
| `--blast_MGEs_minid` | :material-close: | 85 | Coverage (%) threshold to be used when annotating AMR genes with Resfinder, ARGminer and AMRFinderPlus. CARD-RGI is not affected. |
| `--blast_MGEs_mincov` | :material-close: | 85 | Coverage (%) threshold to be used when annotating prophages and mobile elements from PHAST and ICEberg databases |
| `--blast_custom_minid` | :material-close: | 65 | Identity (%) threshold to be used when annotating with user's custom databases |
| `--blast_custom_mincov` | :material-close: | 65 | Coverage (%) threshold to be used when annotating with user's custom databases |

## Merge distance

| <div style="width:200px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--bedtools_merge_distance` | :material-close: | NA | Minimum number of required overlapping bases to merge genes. By default it is not executed. |

## Non-core tools versions

Users can now select the version of the non-core tools Bakta, Unicyler and Flye. These tools now have a parameter which controls which tag, thus version, from quay.io to use.

| Parameter | Default | Description |
| :-------- | :------ | :---------- |
| `--bakta_version`     | 1.7.0--pyhdfd78af_1   | Bakta tool version     |
| `--flye_version`      | 2.9--py39h39abbe0_0   | Flye tool version      |
| `--unicycler_version` | 0.4.8--py38h8162308_3 | Unicycler tool version |

## Max job request options

Set the top limit for requested resources for any single job. If you are running on a smaller system, a pipeline step requesting more resources than are available may cause the Nextflow to stop the run with an error. These options allow you to cap the maximum resources requested by any single job so that the pipeline will run on your system.

!!! note
    
    Note that you can not _increase_ the resources requested by any job using these options. For that you will need your own configuration file. See [the nf-core website](https://nf-co.re/usage/configuration) for details.

| Parameter | Default | Description |
| :-------- | :------ | :---------- |
| `--max_cpus`   | 16     | Maximum number of CPUs that can be requested for any single job   |
| `--max_memory` | 20.GB  | Maximum amount of memory that can be requested for any single job |
| `--max_time`   | 40.h   | Maximum amount of time that can be requested for any single job   |
