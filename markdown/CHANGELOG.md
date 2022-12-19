# bacannot pipeline changelog

The tracking for changes started in v2.1

## v3.2 [19-December-2022]

* Fixes https://github.com/fmalmeida/bacannot/issues/68 reported by @lam-c
    + required a small fix in the scripts used to generated the shiny app
    + thus, new docker image
* Fixes https://github.com/fmalmeida/bacannot/issues/70
    + Adds a `CIRCOS` plot of main annotation results
    + For now it is very minimum and **will** be gradually enhancement in following versions. It contains:
        - forward features
        - reverse features
        - rRNA
        - tRNA
        - AMRFinderPlus and VFDB annotated genes (as labels)
        - PhiSpy annotated prophages
        - GC Skew

> Updates are som many that requires new pipeline feature version, to `3.2` and therefore new docker images.

## v3.1.7 [6-December-2022]

* Fixes https://github.com/fmalmeida/bacannot/issues/75 reported by @sunitj, who highlighted that bakta tool required an update.
* Now, as a consequence, the non-core tools, namely Bakta, Flye and Unicycler that are used only when required in particular cases now have a parameter to control the version of the tool desired.
  * These non-default tools are used with bioconda images from quay.io, now they have parameters `--bakta_version`, `--unicycler_version` and `--flye_version` that allows users to control the version (from quay.io) used.
* PR: https://github.com/fmalmeida/bacannot/pull/76

## v3.1.6 [7-November-2022]

* Fixes https://github.com/fmalmeida/bacannot/issues/71 reported by @fetyj, who highlights that antismash module was failing when outputs for a sample were empty.
  - To solve this issue, the module was updated as such the antismash main results is now optional instead of being required, and 'gff conversion' steps are done only when it's results are not empty.
  - PR: https://github.com/fmalmeida/bacannot/pull/72

## v3.1.5 [17-September-2022]

* Fixes https://github.com/fmalmeida/bacannot/issues/64 reported by @fmalmeida, which highlights that the resfinder annotation were sometimes being duplicated because a single gene had intersection to more then one CDS regions, but, only one being "true".
  - To solve such issue, intersections are now sorted by length, and only the first occurence (per gene) is kept.
* As reported in issue #39, named outputs are now used in modules for better readability.

## v3.1.4 [13-September-2022]

* Fixes https://github.com/fmalmeida/bacannot/issues/62 reported by @rujinlong, where Island-Path tool was failling because it was running on genbank files with no true CDS. This was hapenning because Bakta writes in the comments that the GBK has 0 CDS and, at first, the module was selecting GBK by checking if the CDS string was there. It has now been modified to also work with Bakta.

## v3.1.3 [9-September-2022]

Main changes:

* Hmmer version has been fixed in docker image to allow usage of PGAP HMM on prokka annotation
* Bakta has been added as an option instead of prokka, if desired. Bakta database has to been downloaded manually.
* All the annotation summarisation and reporting modules have been updated to understand when users use Prokka or Bakta.

## v3.1.2

This version contains:

* A fix on the custom blast report file which had an wrong ifelse statement definition
* Updated version of renv docker with fixed report Rmd
* Addition of `checkIfExists` statements when loading file to make sure that user is pointing to files that exists and did not misspelled them
* Fixed a small misspelling on example samplesheet url
* A small fix on phigaro.nf module which was causing the pipeline to fail when its outputs were empty

Thanks to @fetyj for spotting these issues.

## v3.1.1

This version contains a quick fix spotted by @fetyj described in issue [#52](https://github.com/fmalmeida/bacannot/issues/52). Now the path to the argminer backup database is given as full path.

## v3.1

In this version of the pipeline, N things have been added / changed:

1. In order to diminish docker image sizes and to avoid problems with its download hanging and being limited, the databases are not anymore available inside the images. Users must first download the databases with the parameter `--get_dbs`.
    * See https://bacannot.readthedocs.io/en/master/quickstart/#bacannot-databases
2. Since the databases are now downloaded by the pipeline and saved in the user machine, the pipeline now requires the database path as input. E.g. `--bacannot_db ./bacannot_dbs`.
    * See https://bacannot.readthedocs.io/en/master/quickstart/#run-the-pipeline
3. The pipeline now does not have a default profile set, and users must use either `-profile docker` or `-profile singularity`, otherwise it will try to load tools from your system and will surely fail.
4. In previous versions, the pipeline only accepted custom gene databases for additional annotation in nucleotide fasta. Now, the pipeline also accepts protein fastas. Be aware that headers must be properly formatted. Additionally, the pipeline now also accepts a list of `NCBI protein IDs` as additional database with `--ncbi_proteins`. Using it, the pipeline will download the protein sequences and format the databases automatically.
    * See https://bacannot.readthedocs.io/en/master/custom-db/
5. Finally, the pipeline have been reestructured to use and benefit more from the incredible nf-core framework.

## v3.0.1

### additions

In this version of the pipeline, two things have been added:

1. A small dataset with _Haemophilus influenzae_ to provide a `quicktest` profile that runs in 10~20 min.
2. A new github action to test the upcoming changes for each new PR using this available `quicktest`.

### fixes/changes/improvements

These improvements in the code and standardization of namings throughout the pipeline have been made possible thanks to the awesome contributions of [@abhi18av](https://github.com/abhi18av).

1. The config files have been splitout into smaller ones and placed inside a `conf` directory to keep things more organized (PR [#43](https://github.com/fmalmeida/bacannot/pull/43))
2. The names of channels and modules in the pipeline have been changed so they are more standardized and more similar to how it is used for nextflow pipelines in general, so it is easier to read and understand (PR [#45](https://github.com/fmalmeida/bacannot/pull/45))

### comments

> Nothing has changed in terms of how tools are called and used, thus the docker image still the same. In fact, patch/fix releases (x.x.x) will always use the docker from breaking/features release (x.x)

## v3.0

### input configuration

* In order to keeps things the least complex possible and to make the pipeline less confusing, the pipeline has been reconfigured in order to properly use it, in all workflow types (for multiple samples at once or just one) through the samplesheet.
    + Because of that, we removed the possibility to pass the input reads via the command line and now, the files input data files, must always be set inside the samplesheet, even if analysing only one sample.
    + Read more at: https://bacannot.readthedocs.io/en/latest/samplesheet
    + Check the template samplesheet at: https://github.com/fmalmeida/bacannot/blob/master/example_samplesheet.yaml
    + The samplesheet is given with the parameter `--input`
* Due to the implementation above, the folowing parameters are now deprecated, since they are now set inside the YAML file:
    + `--genome`
    + `--sreads_paired`
    + `--sreads_single`
    + `--lreads`
    + `--lreads_type`
    + `--nanopolish_fast5`
    + `--nanopolish_fastq`
* The `--resfinder_species` parameter keeps existing. It now sets a default for all samples in the samplesheet. However, when a sample has another value for that set with the key `resfinder`, the pipeline will use, for that specific sample, the value found inside the samplesheet.

### nomenclature change

* In order to make it simple and natural, two changes ocurred in input/output parameters
    + The `--outdir` parameter is now `--output`
    + The `--in_yaml` parameter is now `--input`
  
### comments

* Since this changes are somewhat major changes, the pipeline main version has changed and it is now in v3.0
    + The docker image is now `fmalmeida/bacannot:v3.0` and `fmalmeida/bacannot:v3.0_renv`

## v2.4.2

Changed how `tag` directives are used inside the pipeline. Now, instead of showing information about the process, it shows which sample is being processed, which is more useful to users.

> Nothing has changed in terms of how tools are called and used, thus the docker image still the same. In fact, patch/fix releases (x.x.x) will always use the docker from breaking/features release (x.x)

## v2.4.1

### hotfix

### hotfix

Super small fix to properly load YAML file when using the pipeline with cloud computing environments such as AWS/S3-bucket:

```bash
# from
parameter_yaml = new FileInputStream(new File(params.in_yaml))
# to
parameter_yaml = file(params.in_yaml).readLines().join("\n")
```

> Nothing has changed in terms of how tools are called and used, thus the docker image still the same. In fact, patch/fix releases (x.x.x) will always use the docker from breaking/features release (x.x)

## v2.4

This release marks the incrementation of the pipeline with two more modules:

1. A new module for the annotation of secondary metabolites with [antiSMASH](https://docs.antismash.secondarymetabolites.org/) has been added.
    + The modules runs only with antiSMASH core modules in order to keep it fast
    + This step can be skipped with `--skip_antismash`
2. A new plugin has been added to the Web (Shiny) App. [SequenceServer](https://sequenceserver.com/) has been implemented inside the app so that users can quickly execute and visualise blast alignments against the samples' genome, genes and proteins.

Also, a small fix was performed in the `run_jbrowse.sh` script. To add the antismash results and to properly check whether digIS results are available or not.

Finally, a small fix in the merge annotation process and SQLparser in the shiny app was also performed to include Resfinder as an option.

Because of that:

+ The `fmalmeida/bacannot:server` has been updated
+ A new image `fmalmeida/bacannot:antismash` is now available for antiSMASH module

## v2.3.2

1. There is a smal fix in the `update_database_image.sh` so it takes only the first two numbers of release tags.
2. Changed the pipeline to now accept the resfinder species panel "Other"
3. Fixed the Resfinder report to understand when only the pointfinder mutations are empty.

## v2.3.1

### Fix in main.nf

Fixed a very small problem that was holding up the execution of flye and unicycler when using the multi-samples workflow with the "samplesheet.yml". The scripts for unicycler and flye under multi-samples workflow was with an "if" statement in the wrong position.

## v2.3

* Dockerfile
    + The dockerfile structure has changed and its size reduced
    + started to track the main tools and databases with dockerfile release tags
    + Github actions wer set to update the main dockerfile containing the databases every month
    + Users can now easily update the databases in the docker image at any time
* New modules
    + A module for identification of top 10 closest Refseq genomes has been added with [RefSeq Masher](https://github.com/phac-nml/refseq_masher)
    + A module for focused detection of insertion sequences has been added with [digIS](https://github.com/janka2012/digIS)
    + A new report has been added to summarise the most general annotation results
* Small fixes
    + The image generated by the vfdb virulence report has been fixed
    + Examples of outputs have been added to the documentation
    + Barrnap is now properly outputting its version

## v2.2

* Automatic reports
    + Some problems with the automatic reports, such as not understanding missing or empty files, and conditionals between processes, have been addressed and fixed.
* VFDB database
    + VFDB database has been set again to use only the core (A) virulence gene dataset. The use of complete (B) dataset were returning lots of spurious hits.
* Resistance annotation
    + The Resfinder software and the ARGminer database have been added to the pipeline in order to augment the diversity of sources for resistance genes annotation.
* Prophage annotation
    + The PhiSpy software has been added to the pipeline to increase the diversity of sources for prophage sequences annotation. This software is complementary to Phigaro.
* KO annotation
    + KeggDecoder has been added to provide a quick and simple visualization of the KEGG Orthology (KO) annotation results.
* Plasmid annotation
    + The Platon software has been added to aid in the prediction of plasmid sequences, together with plasmidfinder.
* Custom annotation
    + The possibility to perform additional annotations using user's custom nucleotide gene databases has been added with the `--custom_db` parameter. The proper configuration of these databases are documented at: https://bacannot.readthedocs.io/en/latest/custom-db
* Multiple genome analysis
    + The possibility to perform the annotation of multiple genomes at once has been added with the `--in_yaml` parameter which takes as input a samplesheet file in YAML format as described in the documentation at: https://bacannot.readthedocs.io/en/latest/samplesheet
* Annotation from raw reads
    + The possibility to annotate genomes from raw reads have been added with the parameters `--sreads_single`, `--sreads_paired` and `--lreads` which takes as input raw sequencing reads from Illumina, Pacbio and ONT in FASTq format and uses it to assemble the genome with Unicycler or Flye (depending on the data availability) prior to the annotation step.
* Bacannot shiny server
    + A simple shiny server has been created and implemented with the `run_server.sh` bash script that loads the shiny app from a docker image that allows the user to quickly interrogate the annotation results via the JBrowse genome browser, the annotation reports and with a BLAST tool implemented in the server that enables users to quickly detect the presence of additional genes/sequences. Take a better look at: https://bacannot.readthedocs.io/en/latest/outputs

## v2.1

This versions have a few additions to the pipeline workflow, they are highlighted and explained below:

* nf-core schema
    + We have added a nextflow parameter schema in json that is compliant with nf-core. This enables that users trigger the graphical interface for configuration and execution of the pipeline via [nf-core launch](https://nf-co.re/launch) utility, also it is possible to track the pipeline execution with [nextflow tower](https://tower.nf/).

```bash
# It is triggered as
nf-core launch fmalmeida/bacannot
```

* plasmidfinder
    + _In silico_ plasmid detection has been added with Plasmidfinder. In order to not execute this parameter, users will need to use `--not_run_plasmid_search`. Otherwise, its thresholds can be configured with `--plasmids_minid` and `--plasmids_mincov`.
* VFDB database
    + VFDB database has been changed from the core (A) dataset to the complete (B) dataset. This dataset is first clustered with cd-hit using a 90% identity threshold and this new database file (after cd-hit) is used for the virulence gene annotation.
