process VICTORS_DB {
    publishDir "${params.output}/victors_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_ultralow' ]

    output:
    file("*")

    script:
    """
    # download victors database (aa)
    ( 
        wget --tries=10 -O victors_original.fasta "http://www.phidias.us/victors/downloads/gen_downloads_protein.php"
        grep -v "^[^>M]" victors_original.fasta > victors_prot.fasta && \\
            rm victors_original.fasta && \\
            awk -v db=victors '/>/{ split(\$0,a,"|"); split(a[5],gene," \\\\["); all=\$0; \$0=">" db "~~~" gene[1] "~~~" a[4] "~~~" "Victors_" a[2] "~~~" all }1' victors_prot.fasta | \\
            sed -e 's/ >/ /g' -e 's/~ /~/g' | \\
            awk -F "~~~" ' { if (\$0 ~ />/) { gsub(" ", "_", \$2); gsub(" ", "_", \$5); print \$1 "~~~" \$2 "~~~" \$3 "~~~" \$4 "~~~" \$5 } else { print \$0 }}' | \\
            sed -e 's/~~~>/~~~/g' -e 's/|_/|/g' > sequences && \\
            diamond makedb --in sequences -d diamond && \\
            makeblastdb -in sequences -title 'victors' -dbtype prot -logfile /dev/null && \\
            rm victors_prot.fasta
    ) || 
    ( 
        cat /work/victors.fasta > sequences && \\
            diamond makedb --in sequences -d diamond && \\
            makeblastdb -in sequences -title 'victors' -dbtype prot -logfile /dev/null && \\
            rm victors_prot.fasta
    )
    """
}
