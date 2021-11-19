process PLATON_DB {
    publishDir "${params.output}/platon_db", mode: 'copy', overwrite: "$params.force_update"
   
    output:
    file("*")

    script:
    """   
    # download platon database
    wget -O db.tar.gz "https://zenodo.org/record/4066768/files/db.tar.gz"
   """
}
