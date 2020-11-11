# bacannot pipeline changelog

The tracking for changes started in v2.1

## v2.1

This versions have a few additions to the pipeline workflow, they are highlighted and explained below:

### nf-core schema

We have added a nextflow parameter schema in json that is compliant with nf-core. This enables that users trigger the graphical interface for configuration and execution of the pipeline via [nf-core launch](https://nf-co.re/launch) utility, also it is possible to track the pipeline execution with [nextflow tower](https://tower.nf/).

```bash
# It is triggered as
nf-core launch fmalmeida/bacannot
```

Checkout the paremeters `--use_tower` and `--tower_token` to activate pipeline execution in nextflow tower.

### plasmidfinder

_In silico_ plasmid detection has been added with Plasmidfinder. In order to not execute this parameter, users will need to use `--not_run_plasmid_search`. Otherwise, its thresholds can be configured with `--plasmids_minid` and `--plasmids_mincov`.

### VFDB database

VFDB database has been changed from the core (A) dataset to the complete (B) dataset. This dataset is first clustered with cd-hit using a 90% identity threshold and this new database file (after cd-hit) is used for the virulence gene annotation.
