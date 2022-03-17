# Selecting between profiles

## What are profiles?

Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string. They are a set of "sensible defaults" for the resource requirements of each of the steps in the workflow, that can be enabled with the command line flag `-profile`. You can learn more about nextflow profiles at:

+ <https://nf-co.re/usage/configuration#basic-configuration-profiles>
+ <https://www.nextflow.io/docs/latest/config.html#config-profiles>

## Available profiles

The pipeline have "standard profiles" set to run the workflows with either **conda**, **docker** or **singularity** using the [local executor](https://www.nextflow.io/docs/latest/executor.html), which is nextflow's default and basically runs the pipeline processes in the computer where Nextflow is launched.

If you need to run the pipeline using another executor such as sge, lsf, slurm, etc. you can take a look at [nextflow's manual page](https://www.nextflow.io/docs/latest/executor.html) to proper configure one in a new custom profile set in your personal copy of [MpGAP config file](https://github.com/fmalmeida/phylogram/blob/master/nextflow.config) and take advantage that nextflow allows multiple profiles to be used at once, e.g. `-profile conda,sge`.

Please, also check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

!!! note

    If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. **This is not recommended**.

### Note on sigularity

If you are using `singularity` and are persistently observing issues downloading Singularity images directly due to timeout or network issues, try downloading it first. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.

```bash
# for singularity
# remember to properly set NXF_SINGULARITY_LIBRARYDIR
# read more at https://www.nextflow.io/docs/latest/singularity.html#singularity-docker-hub
export NXF_SINGULARITY_LIBRARYDIR=MY_SINGULARITY_IMAGES    # your singularity storage dir
export NXF_SINGULARITY_CACHEDIR=MY_SINGULARITY_CACHE       # your singularity cache dir
singularity pull --dir $NXF_SINGULARITY_LIBRARYDIR fmalmeida-dev.img docker://fmalmeida/dev

# run
nextflow run fmalmeida/phylogram -profile singularity [OPTIONS]
```

### Note on conda

If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.