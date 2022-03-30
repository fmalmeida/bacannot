process ANTISMASH_DB {
    publishDir "${params.output}/antismash_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_ultralow' ]
   
    output:
    file("*")

    script:
    """
    # download antismash database
    export PATH=/opt/conda/envs/antismash/bin:\$PATH
    download-antismash-databases --database-dir \$(pwd)
    """
}
