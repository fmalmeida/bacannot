process MLST_DB {
    publishDir "${params.output}/mlst_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_ultralow' ]

    output:
    file("*")

    script:
    """
    # download mlst database
    curl https://pubmlst.org/static/data/dbases.xml > dbases.xml
    """
}
