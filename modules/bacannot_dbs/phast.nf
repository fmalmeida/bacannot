process PHAST_DB {
    publishDir "${params.output}/phast_db", mode: 'copy', overwrite: "$params.force_update"
    label 'db_download'
   
    output:
    file("*")

    script:
    """
    # download phast database (aa)
    wget -O phast_prot.fasta http://phaster.ca/downloads/prophage_virus.db && \\
        awk -v db=phast '/>/{ split(\$0,a,"|"); split(a[5],gene," \\\\["); all=\$0; \$0=">" db "~~~" gene[1] "~~~" a[4]"~~~" "PHAST_" a[2] " " all }1' phast_prot.fasta | \\
        sed -e 's/ >/ /g' -e 's/~ /~/g' | \\
        awk -F "~~~" ' { if (\$0 ~ />/) { gsub(" ", "_", \$2); print \$1 "~~~" \$2 "~~~" \$3 "~~~" \$4 } else { print \$0 }}' | \\
        awk -F "~~~" ' { if (\$0 ~ />/) { gsub("-", "_", \$2); print \$1 "~~~" \$2 "~~~" \$3 "~~~" \$4 } else { print \$0 }}' > sequences && \\
        rm phast_prot.fasta && \\
        diamond makedb --in sequences -d diamond && \\
        makeblastdb -in sequences -title 'phast' -dbtype prot -logfile /dev/null
    """
}