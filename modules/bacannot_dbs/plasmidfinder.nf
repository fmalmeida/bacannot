process PLASMIDFINDER_DB {
    publishDir "${params.output}", mode: 'copy', overwrite: "$params.force_update"
    label 'db_download'
   
    output:
    file("*")

    script:
    """   
    # download plasmidfinder database
    git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git
    """
}