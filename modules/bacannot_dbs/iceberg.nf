process ICEBERG_DB {
    publishDir "${params.output}/iceberg_db", mode: 'copy', overwrite: "$params.force_update"
    label 'db_download'
   
    output:
    file("*")

    script:
    """
    # download iceberg database (nt)
    wget https://bioinfo-mml.sjtu.edu.cn/ICEberg2/download/ICE_seq_experimental.fas && \\
        awk -v db=ICEberg '/>/{ split(\$0,a,"|"); all=\$0; \$0=">" db "~~~" "ICE_" a[2] "~~~" a[5] "~~~" a[3] " " all }1' ICE_seq_experimental.fas | \\
        sed -e 's/ >/ /g' > sequences && \\
        rm ICE_seq_experimental.fas && \\
        makeblastdb -in sequences -title 'iceberg' -dbtype nucl -logfile /dev/null

    # download iceberg database (aa)
    wget https://bioinfo-mml.sjtu.edu.cn/ICEberg2/download/ICE_aa_experimental.fas && \\
        awk -v db=ICEberg \\
        '/>/{ split(\$0,col," "); split(col[1],a,"[|]"); split(col[2],b,"[|]"); split(\$0,c,"[|]"); all=\$0; \$0=">" db "~~~" "ICE_" a[2] "~~~" b[4] "~~~" c[6] " " all }1' \\
        ICE_aa_experimental.fas | sed -e 's/ >/ /g' | awk -F '\\\\]' \\
        '{ if (\$0 ~ />/) { gsub(" ","_",\$1); gsub("_\\\\[","_",\$1); gsub("~_","~",\$1); print \$1,\$2 "]" } else { print \$0 }}' > proteins && \\
        diamond makedb --in proteins -d diamond && \\
        makeblastdb -in proteins -title 'iceberg' -dbtype prot -logfile /dev/null && \\
        rm ICE_aa_experimental.fas
    """
}
