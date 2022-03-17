# Quickstart

To provide a meaningful quickstart we will use the 31 TNL NLR sequences from [_Actinidia eriantha_](https://plants.ensembl.org/Arabidopsis_thaliana/Info/Index) public data hosted in [ANNA database](https://biobigdata.nju.edu.cn/ANNA/) to generate a quick plot.

## Required inputs

To run the pipeline, we basically need a FASTA of proteins that we want to run the phylogenetic analysis with (`--fasta`) and a Pfam HMM domain which we want to use as the core for aligning and reconstructing phylogeny (`--hmm`).

## Downloading the inputs

### Input sequences

```bash
# get the NLR proteins
wget https://github.com/fmalmeida/test_datasets/raw/main/phylogram_testing/ANNA_ericales_TNL.fasta
```

!!! tip "Input sequences headers"

    The input sequences may contain metadata about the annotation and source of features. For example, the sequence headers from ANNA look like this:

    ```bash
    >DTZ79_22g05730 TNL Actinidia_eriantha
    >DTZ79_29g05250 TNL Actinidia_eriantha 
    [...]
    ```

    This metadata from the headers are used to generate useful annotation files for [Itol](https://itol.embl.de/). Please read more about giving metadata to input sequences [here](seqmetadata.md#).

### Pfam database

The pipeline is capable of generating annotation files for visualizing the phylogeny with domains annotation in Itol. Howeverm for that, we need to provide the Pfam-A.hmm file to the pipeline. Thus, we need to download it:

```bash
# get Pfam-A.hmm
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
```

## Selecting phylogeny references

Besides the required inputs, we can also select NLR reference sequences from [RefPlantNLR collection](https://zenodo.org/record/3936022#.Yij6UsvMJMg) to be used in phylogeny with `--refplants_species`, if desired. A list of options can be found [here](https://github.com/fmalmeida/nlrpipe/blob/main/assets/refplants/species.txt).

!!! tip "Selecting refplants species"

    RefPlants species can be selected in species or genus level. For example, both uses are correct:

    1. `--refplants_species` "Capsicum chacoense"
        * gets only _Capsicum chacoense_ sequences
    2. `--refplants_species` "Capsicum"
        * gets all _Capsicum_ sequences

In our example, to make things simple, let's use the sequences from _Arabidopsis_.

## Running the pipeline

To run our quickstart we would need something like this:

``` { .bash .annotate hl_lines="6 8" }
# running quickstart
nextflow run fmalmeida/phylogram \
    -r main -profile docker \
    --prefix "actinidia_eriantha" \
    --fasta "ANNA_ericales_TNL.fasta" \
    --hmm "PF00931" \
    --refplants_species "Arabidopsis" \
    --bootstraps 2 \
    --pfam_db Pfam-A.hmm
```

* `--hmm`
    * Since our sequences are NLR, let's use the NB-ARC as the core for alignments. Users can select and HMM either giving path to the hmm file in his computer or giving the Pfam ID. If given the Pfam ID, the pipeline will download it.
* `--bootstraps`
    * Let's have IQTree to generate a quick tree by having only two bootstraps iterations. Remember, to have a consensus tree, the minimum value accepted for this parameter is 2.

## Outputs

After a successfull run you will have an output directory like this:

```bash
results
├── CLEAN_HEADER
├── CLIPKIT
├── concatenated_input.faa
├── DOMAIN_ALIGNMENT
├── FASTTREE
├── GET_REFPLANTS
├── IQTREE
├── ITOL_ANNOTATION
├── pipeline_info
├── TCOFFEE_TCS
├── TREESHRINK
└── TRIMAL
```

These are all the main/core files that have been generated during the analysis. The main results are under `FASTTREE` and `IQTREE` or `RAXML` directory which contain the phylogeneyic trees generated. Inside `TREESHRINK` are the pruned (long branches removed) phylogenetic trees.

Finally, inside `ITOL_ANNOTATION` are the annotation files (itol_annotation.*.txt) that can be loaded in Itol together with the phylogenetic tree for visualizing an annotated tree. The `FASTTREE` tree from this quickstart is available in Itol for visualization: <https://itol.embl.de/tree/189635130325201647471973>

> In the Itol annotation: label colors are different sequence sources, for example, different organisms, and the label background color are the different NLR classes.

<img src="../ericales_example_tree.svg" width="90%">