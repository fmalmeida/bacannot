# Installation

## Dependencies

The pipeline require only a UNIX system, [Nextflow](https://www.nextflow.io/docs/latest/index.html#) and either [Docker](https://www.docker.com/), [Singularity](https://sylabs.io/docs/) or [Conda](https://docs.conda.io/en/latest/). Please, for installing these tools refer to their manual.

## Downloading the pipeline

You can easily get a copy of the pipeline with:

```bash
# nextflow pull
nextflow pull fmalmeida/phylogram
```

## Testing your installation

After that, you can run the pipeline with a testing dataset by selecting one of the available profiles: 

1. Docker
    * `nextflow run fmalmeida/mpgap -profile docker,test`
2. Singularity
    * `nextflow run fmalmeida/mpgap -profile singularity,test`
3. Conda
    * `nextflow run fmalmeida/mpgap -profile conda,test`

!!! note "About NF profiles"

    Please read more about how to [proper select NF profiles](https://github.com/fmalmeida/phylogram#selecting-between-profiles) to better understand it.