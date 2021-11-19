process CREATE_DBS {
    publishDir "${params.output}", mode: 'copy', overwrite: "$params.force_update"
   
    output:
    file("*")

    script:
    """
    # download iceberg database (nt)
    mkdir iceberg_db && \\
        wget https://bioinfo-mml.sjtu.edu.cn/ICEberg2/download/ICE_seq_experimental.fas && \\
        awk -v db=ICEberg '/>/{ split(\$0,a,"|"); all=\$0; \$0=">" db "~~~" "ICE_" a[2] "~~~" a[5] "~~~" a[3] " " all }1' ICE_seq_experimental.fas | \\
        sed -e 's/ >/ /g' > iceberg_db/sequences && \\
        rm ICE_seq_experimental.fas && \\
        makeblastdb -in iceberg_db/sequences -title 'iceberg' -dbtype nucl -logfile /dev/null

    # download iceberg database (aa)
    wget https://bioinfo-mml.sjtu.edu.cn/ICEberg2/download/ICE_aa_experimental.fas && \\
        awk -v db=ICEberg \\
        '/>/{ split(\$0,col," "); split(col[1],a,"[|]"); split(col[2],b,"[|]"); split(\$0,c,"[|]"); all=\$0; \$0=">" db "~~~" "ICE_" a[2] "~~~" b[4] "~~~" c[6] " " all }1' \\
        ICE_aa_experimental.fas | sed -e 's/ >/ /g' | awk -F '\\]' \\
        '{ if (\$0 ~ />/) { gsub(" ","_",\$1); gsub("_\\[","_",\$1); gsub("~_","~",\$1); print \$1,\$2 "]" } else { print \$0 }}' > iceberg_db/proteins && \\
        diamond makedb --in iceberg_db/proteins -d iceberg_db/diamond && \\
        makeblastdb -in iceberg_db/proteins -title 'iceberg' -dbtype prot -logfile /dev/null && \\
        rm ICE_aa_experimental.fas

    # download victors database (aa)
    mkdir victors_db && \
        wget -O victors_original.fasta "http://www.phidias.us/victors/downloads/gen_downloads_protein.php" && \\
        grep -v "^[^>M]" victors_original.fasta > victors_prot.fasta && \\
        rm victors_original.fasta && \\
        awk -v db=victors '/>/{ split(\$0,a,"|"); split(a[5],gene," \\["); all=\$0; \$0=">" db "~~~" gene[1] "~~~" a[4] "~~~" "Victors_" a[2] " " all }1' victors_prot.fasta | \\
        sed -e 's/ >/ /g' -e 's/~ /~/g' | \\
        awk -F "~~~" ' { if (\$0 ~ />/) { gsub(" ", "_", \$2); print \$1 "~~~" \$2 "~~~" \$3 "~~~" \$4 } else { print \$0 }}' > victors_db/sequences && \\
        diamond makedb --in victors_db/sequences -d victors_db/diamond && \\
        makeblastdb -in victors_db/sequences -title 'victors' -dbtype prot -logfile /dev/null && \\
        rm victors_prot.fasta

    # download phast database (aa)
    mkdir phast_db && \
        wget -O phast_prot.fasta http://phaster.ca/downloads/prophage_virus.db && \\
        awk -v db=phast '/>/{ split(\$0,a,"|"); split(a[5],gene," \\["); all=\$0; \$0=">" db "~~~" gene[1] "~~~" a[4]"~~~" "PHAST_" a[2] " " all }1' phast_prot.fasta | \\
        sed -e 's/ >/ /g' -e 's/~ /~/g' | \\
        awk -F "~~~" ' { if (\$0 ~ />/) { gsub(" ", "_", \$2); print \$1 "~~~" \$2 "~~~" \$3 "~~~" \$4 } else { print \$0 }}' | \\
        awk -F "~~~" ' { if (\$0 ~ />/) { gsub("-", "_", \$2); print \$1 "~~~" \$2 "~~~" \$3 "~~~" \$4 } else { print \$0 }}' > phast_db/sequences && \\
        rm phast_prot.fasta && \\
        diamond makedb --in phast_db/sequences -d phast_db/diamond && \\
        makeblastdb -in phast_db/sequences -title 'phast' -dbtype prot -logfile /dev/null

    # download argminer database (aa)
    ## argminer server has a lot of problems
    mkdir argminer_db && \\
        wget -t 1 https://github.com/fmalmeida/bacannot/raw/master/docker/argminer_bkp/argminer.fasta && \\
        awk -v db=ARGMiner '/>/{ split(\$0,a,"|"); \$0=">" db "~~~" a[3] "~~~" a[1] "~~~" a[2] " " a[4] }1' argminer.fasta | \\
        sed -e 's/~>/~/g' -e 's/gi:.*:ref://g' -e 's/gi:.*:gb://g' -e 's/gi:.*:emb://g' -e 's/:~/~/g' > argminer_db/sequences && \\
        rm argminer.fasta && \\
        makeblastdb -in argminer_db/sequences -title 'argminer' -dbtype prot -logfile /dev/null && \\
        diamond makedb --in argminer_db/sequences -d argminer_db/diamond
    """
}
