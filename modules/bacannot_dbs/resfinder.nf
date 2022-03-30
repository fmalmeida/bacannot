process RESFINDER_DB {
    publishDir "${params.output}/resfinder_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_ultralow' ]
   
    output:
    file("*")

    script:
    """
    # download resfinder databases
    git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git db_resfinder
    rm -r db_resfinder/.git
    git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git db_pointfinder
    rm -r db_pointfinder/.git
    """
}
