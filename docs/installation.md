# Installation

## Dependencies

The pipeline require only a UNIX system, [Nextflow](https://www.nextflow.io/docs/latest/index.html#) and either [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/docs/). Please, for installing these tools refer to their manual.

## Downloading the pipeline

You can easily get a copy of the pipeline with:

```bash
# nextflow pull
nextflow pull fmalmeida/bacannot
```

!!! warning
    
    The pipeline requires a UNIX system, therefore, Windows users may successfully use this pipeline via the [Linux subsystem for window](https://docs.microsoft.com/pt-br/windows/wsl/install-win10). Nextflow team has made available a [nice tutorial](https://www.nextflow.io/blog.html) about this issue.

## Downloading docker images

> The pipeline uses both custom and public images.
> All images can be downloaded on the fly, automatically by nextflow, and this is the recommended way to do it.

If you want to download it yourself, you can find all the images used in the pipeline described in the file [docker.config](https://github.com/fmalmeida/bacannot/blob/master/conf/docker.config) (for docker) and [singularity.config](https://github.com/fmalmeida/bacannot/blob/master/conf/singularity.config) (for singularity).

The images are defined like the following:

```bash
...
withLabel: 'db_download|db_tools|misc' {
    container = 'fmalmeida/bacannot@sha256:0648797837cd8e11b6abd40745cafc0db81647953921ec54ce0ceef9ecef6450'
}
...
```

And could be downloaded like this:

```bash
docker pull fmalmeida/bacannot@sha256:0648797837cd8e11b6abd40745cafc0db81647953921ec54ce0ceef9ecef6450
```

> You would need to do it for each image.

!!! info "If using singularity"

    **Docker and singularity images are downloaded on the fly**. Be sure to properly set `NXF_SINGULARITY_LIBRARYDIR` env variable to a writable directory if using Singularity. This will make that the downloaded images are reusable through different executions. Read more at: https://www.nextflow.io/docs/latest/singularity.html#singularity-docker-hub

    For example, to download the images for docker you may:

    ```bash
    # apply this command to each image
    # just change the "/" and ":" for "-".
    # E.g. Image fmalmeida/bacannot:v3.3_misc becomes fmalmeida-bacannot-v3.3_misc.img
    singularity pull --dir $NXF_SINGULARITY_LIBRARYDIR fmalmeida-bacannot-v3.3_misc.img docker://fmalmeida/bacannot:v3.3_misc
    ```

## Bacannot databases

Bacannot databases are not inside the docker images anymore to avoid huge images and problems with connections and limit rates with dockerhub.

### Pre-formatted

Users can directly download pre-formatted databases from Zenodo: https://doi.org/10.5281/zenodo.7615811

Useful for standardization and also overcoming known issues that may arise when formatting databases with `singularity` profile.

A module to download the latest pre-formatted database has also been made available:

```bash
# Download pipeline pre-built databases
nextflow run fmalmeida/bacannot \
    --get_zenodo_db \
    --output ./ \
    -profile <docker/singularity>
```

### I want to generate a new formatted database

```{bash .annotate hl_lines="5"}
# Download pipeline databases
nextflow run fmalmeida/bacannot \
    --get_dbs \
    --output bacannot_dbs \
    -profile <docker/singularity>
```

## Testing your installation

After that, you can run the pipeline with a testing dataset by selecting one of the available profiles: 

1. Docker
    * `nextflow run fmalmeida/mpgap -profile docker,test` --bacannot_db ./bacannot_dbs
2. Singularity
    * `nextflow run fmalmeida/mpgap -profile singularity,test` --bacannot_db ./bacannot_dbs

!!! note "About NF profiles"

    Please read more about how to [proper select NF profiles](profiles.md#) to better understand it.
