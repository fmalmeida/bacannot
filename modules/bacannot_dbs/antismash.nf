process ANTISMASH_DB {
    publishDir "${params.output}/antismash_db", mode: 'copy', overwrite: "$params.force_update"
    label 'db_download'
   
    output:
    file("*")

    script:
    """
    # download antismash database
    set +eu ; source activate antismash
    download-antismash-databases --database-dir \$(pwd)
    conda deactivate
    """
}
