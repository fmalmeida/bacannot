# Runnig earlier releases

By default, nextflow will execute the code that is available in master/main branch of a repository. This is the recommended action because the main branch will always hold the most up-to-date code with bug fixes from other versions.

However, if for some reason you need to execute a different version of a pipeline, you can use the parameter `-r branch/tag`. For example:

```bash
# to run the code in tag/branch v2.0
nextflow run fmalmeida/bacannot -r v2.0 --help
```

However, in the first releases of the pipeline, I was still new to the creation and management of github releases. Sometimes I found a bug, fixed and commited to main without triggering a small release to contain the fixes. Also, was not using frozen docker images to specific version releases.

Therefore, if you simply run `nextflow run fmalmeida/bacannot -r v2.0` the pipeline will probably fail because the code in the first release tags (v2.0, v2.1, v2.2) were not keeping track of small bug fixes. This is my fault.

So, for these releases, I have created 3 branches to store latest codes (with bug fixes) for each one of this releases. They are called `{v2.0,v2.1,v2.2}_compatibility`.

So, to run these specific versions, users must use:

```bash
# to run the code in tag/branch v2.1
nextflow run fmalmeida/bacannot -r v2.1_compatibility --help
```

## Note

I have now understood how to best use and manage the github releases. Thus, for further releases, small fixes will be tracked with small releases (such as v2.3.1, v2.3.2, etc.) and, each release will have a frozen docker image.

Thus, for releases from v2.3 and beyond, the release tags could be used smoothly as held by nextflow as: `nextflow run fmalmeida/bacannot -r {release tag}`.
