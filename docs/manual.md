# Manual

## Input/output options

| <div style="width:50px">Parameter</div> | Required | Default | Description |
| :-------------------------------------- | :------- | :------ | :---------- |
| `--fasta` | :material-check: | NA | Input FASTA (one or many) of peptides for analysis |
| `--prefix` | :material-close: | input | Prefix for output files |
| `--output` | :material-close: | results | Output directory where the results will be saved |

!!! note "About input fasta"

    If more than one is given. they will be concatenated and analyzed together. Your input may have metadata, please read more [here](seqmetadata.md#).

    The "^" character is not allowed because of Itol and is automatically replaced by "-". The input sequences with fixed headers will be available in the root output directory (`concatenated_input.fasta`).

    Also, even though it does not cause problems, beware that by default Itol does not show "_" characteres, so you may want to use "-" instead.

## Reference/Model options

| <div style="width:160px">Parameter</div> | Required | Default | Description |
| :--------------------------------------- | :------- | :------ | :---------- |
| `--hmm` | :material-check: | PF00931 | Path to desired HMM model file or Pfam ID. If Pfam ID is given, the pipeline will download the model. |
| `--aln_min_len` | :material-check: | 20 | Min. length of sequence alignments against HMM domain to keep |
| `--refplants_species` | :material-close: | NA | Select species to load reference RefPlants sequences. A list of options can be found [here](https://github.com/fmalmeida/nlrpipe/blob/main/assets/refplants/species.txt). |
| `--pfam_db` | :material-close: | NA | Path to [Pfam-A.hmm](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz) file. This file is required to run the Itol annotation module. If not given, Itol files will not be generated. |
| `--domains_metadata` | :material-close: | [NB-ARC template](https://github.com/fmalmeida/phylogram/blob/main/assets/NB-ARC_Domains_Itol_Guide.tsv) | This file is a three column metadata of the main pfam domains to be parsed for Itol legends |

!!! note "About domains metadata"

    By default, the Itol Annotation step checks for all good quality Pfam Domains annotation. But, shapes and colors are defined randomly and, since sometimes many different domains may appear, they are not drawn into Itol legend by default. And this is the purpose of this parameter!
    
    The TSV has three columns (**Domain_Name**    **Itol_Shape**    **Itol_Color**). This is setup which domains for the annotation should be considered the core of the analyses. If present in this file, instead of generating randomly, the pipeline will use the shapes and colors provided in the file and, generate a legend for Itol using the entries provided. Here it is [the template](https://github.com/fmalmeida/phylogram/blob/main/assets/NB-ARC_Domains_Itol_Guide.tsv) used by the pipeline it none is provided.

!!! warning "About Min. Alignment Length"

    Hmmsearch sequence alignments are filtered based on minimum alignment length with easel tools using the value set with `--aln_min_len`. Sometimes your sequence may not appear in the final phylogeny because it was filtered out from MSA due to the minimum alignment length allowed. If this happens, you may try lowering this parameter. In theory, if set to 1, anything that aligns will be allowed.

## Phylogeny options

| <div style="width:160px">Parameter</div> | Default | Description |
| :--------------------------------------- | :------ | :---------- |
| `--use_raxml` | false | Use RAXML instead of IQTREE to compute phylogeny. If using RAXML, the pipeline will use [jmodeltest](https://github.com/ddarriba/jmodeltest2) to select best model |
| `--bootstraps` | 100 | How many bootstraps do you desire the phylogeny software to run? |
| `--use_fastbootstrap` | false | Use [iqtree's fast bootstrap mode](http://www.iqtree.org/doc/Tutorial#assessing-branch-supports-with-ultrafast-bootstrap-approximation)? Option only useful for IQTREE |
| `--use_modeltest` | false | Use [jmodeltest instead of iqtree's modelfinder](http://www.iqtree.org/doc/Tutorial#choosing-the-right-substitution-model)? Option only useful for IQTREE |

!!! note "About fastbootstrap"

    IQTree does not allow the use of fastbootstrap for less than 1000 bootstrap iterations. In addition, we know that user's may not want to always use the fastbootstrap algorithm.

    Therefore, in order to avoid having different behaviours depending on the amount of bootstraps we set the pipeline to do not use IQTree's fastbootstrap algorithm by default.

    Users must explicitely ask for it with `--use_fastboostrap`. Beware that if your `--bootstraps` is lower than 1000, when using `--use_fastbootstrap` it will be automatically changed to 1000. 

## Workflow modules management

| Parameter | Default | Description |
| :-------- | :------ | :---------- |
| `--skip_clipkit` | false | Deactivate Clipkit module     |
| `--skip_trimal`  | false | Deactivate Trimal module      |
| `--skip_tcoffee` | false | Deactivate Tcoffee TCS module |

!!! note "About tcoffee TCS"

    [Transitive Consistency Score (TCS)](https://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#transitive-consistency-score-tcs)

    TCS is an alignment evaluation score that makes it possible to identify the most correct positions in an MSA. It has been shown that these positions are the most likely to be structuraly correct and also the most informative when estimating phylogenetic trees. The TCS evaluation and filtering procedure is implemented in the T-Coffee package and can be used to evaluate and filter any third party MSA (including T-Coffee MSA of course!).

## Max job request options

Set the top limit for requested resources for any single job. If you are running on a smaller system, a pipeline step requesting more resources than are available may cause the Nextflow to stop the run with an error. These options allow you to cap the maximum resources requested by any single job so that the pipeline will run on your system.

!!! note
    
    Note that you can not _increase_ the resources requested by any job using these options. For that you will need your own configuration file. See [the nf-core website](https://nf-co.re/usage/configuration) for details.

| Parameter | Default | Description |
| :-------- | :------ | :---------- |
| `--max_cpus`   | 16     | Maximum number of CPUs that can be requested for any single job   |
| `--max_memory` | 128.GB | Maximum amount of memory that can be requested for any single job |
| `--max_time`   | 240.h  | Maximum amount of time that can be requested for any single job   |
