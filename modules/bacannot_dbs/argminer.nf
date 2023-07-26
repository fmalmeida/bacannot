process ARGMINER_DB {
    publishDir "${params.output}/argminer_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_ultralow' ]

    output:
    file("*")

    script:
    """
    # download argminer database (aa)
    ## argminer server has a lot of problems
    ( 
        wget -t 1 http://bench.cs.vt.edu/ftp/argminer/release/ARGminer-v1.1.1.A.fasta && \\
		awk -v db=ARGMiner '/>/{ split(\$0,a,"|"); \$0=">" db "~~~" a[3] "~~~" a[1] "~~~" a[2] "~~~" a[4] }1' ARGminer-v1.1.1.A.fasta | \\
		sed -e 's/~>/~/g' -e 's/gi:.*:ref://g' -e 's/gi:.*:gb://g' -e 's/gi:.*:emb://g' -e 's/:~/~/g' -e 's/:_/_/g' -e 's/ /_/g' > sequences && \\
		rm ARGminer-v1.1.1.A.fasta && \\
		makeblastdb -in sequences -title 'argminer' -dbtype prot -logfile /dev/null && \\
		diamond makedb --in sequences -d diamond
    ) || 
    ( 
        cat /work/argminer.fasta | \\
        sed -e 's/~>/~/g' -e 's/gi:.*:ref://g' -e 's/gi:.*:gb://g' -e 's/gi:.*:emb://g' -e 's/:~/~/g' -e 's/:_/_/g' -e 's/ /_/g' > sequences && \\
		makeblastdb -in sequences -title 'argminer' -dbtype prot -logfile /dev/null && \\
		diamond makedb --in sequences -d diamond 
    )
    """
}
