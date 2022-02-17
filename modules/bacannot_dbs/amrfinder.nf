process AMRFINDER_DB {
    publishDir "${params.output}/amrfinder_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_ultralow' ]
   
    output:
    file("*")

    script:
    """
    # download amrfinderplus database
    amrfinder_update -d \$(pwd)
    """
}
