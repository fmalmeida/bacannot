process RESFINDER_DB {
    publishDir "${params.output}/resfinder_db", mode: 'copy', overwrite: "$params.force_update"
   
    output:
    file("*")

    script:
    """
    # download resfinder databases
    git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git db_resfinder
    git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git db_pointfinder
    """
}
