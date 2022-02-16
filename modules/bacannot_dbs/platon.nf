process PLATON_DB {
    publishDir "${params.output}/platon_db", mode: 'copy', overwrite: "$params.force_update"
    label 'db_download'
   
    output:
    file("*")

    script:
    """   
    # download platon database
    wget -O db.tar.gz "https://zenodo.org/record/4066768/files/db.tar.gz"
    tar zxvf db.tar.gz
    mv db/* .
    rm -rf db db.tar.gz
    """
}
