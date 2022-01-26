process ARGMINER_DB {
    publishDir "${params.output}/argminer_db", mode: 'copy', overwrite: "$params.force_update"
    label 'db_download'
   
    output:
    file("*")

    script:
    """
    # download argminer database (aa)
    ## argminer server has a lot of problems
    wget https://github.com/fmalmeida/bacannot/raw/master/docker/argminer_bkp/argminer.fasta -O sequences && \\
        makeblastdb -in sequences -title 'argminer' -dbtype prot -logfile /dev/null && \\
        diamond makedb --in sequences -d diamond
    """
}
