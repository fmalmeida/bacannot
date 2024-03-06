process ANTISMASH_DB {
    publishDir "${params.output}/antismash_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_ultralow' ]

    output:
    file("*")

    script:
    def antismash_version='6.1.1'

    if (params.running_engine == 'singularity')
    """
    mkdir local-install
    export PYTHONUSERBASE=./local-install
    export PATH=/opt/conda/envs/antismash/bin:\$PATH

    # install locally so it can download dbs
    # singularity has many read-write permissions for this tool
    wget https://dl.secondarymetabolites.org/releases/${antismash_version}/antismash-${antismash_version}.tar.gz
    tar zxvf antismash-${antismash_version}.tar.gz
    python -m pip install --user ./antismash-${antismash_version}
    export PYTHONPATH=\$(realpath \$( find ./local-install -name 'site-packages' ))

    # now download it
    # download antismash database
    ./local-install/bin/download-antismash-databases --database-dir ./

    # delete it
    rm -rf ./local-install ./antismash-${antismash_version}
    """

    else
    """
    # download antismash database
    export PATH=/opt/conda/envs/antismash/bin:\$PATH
    download-antismash-databases --database-dir ./
    """
}
