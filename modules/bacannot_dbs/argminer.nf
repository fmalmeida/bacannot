process ARGMINER_DB {
    publishDir "${params.output}/argminer_db", mode: 'copy', overwrite: "$params.force_update"
   
    output:
    file("*")

    script:
    """
    # download argminer database (aa)
    ## argminer server has a lot of problems
    wget https://github.com/fmalmeida/bacannot/raw/master/docker/argminer_bkp/argminer.fasta && \\
        awk -v db=ARGMiner '/>/{ split(\$0,a,"|"); \$0=">" db "~~~" a[3] "~~~" a[1] "~~~" a[2] " " a[4] }1' argminer.fasta | \\
        sed -e 's/~>/~/g' -e 's/gi:.*:ref://g' -e 's/gi:.*:gb://g' -e 's/gi:.*:emb://g' -e 's/:~/~/g' > sequences && \\
        rm argminer.fasta && \\
        makeblastdb -in sequences -title 'argminer' -dbtype prot -logfile /dev/null && \\
        diamond makedb --in sequences -d argminer_db/diamond
    """
}
