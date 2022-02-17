process VFDB_DB {
    publishDir "${params.output}/vfdb_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_ultralow' ]
   
    output:
    file("*")

    script:
    """
    # download vfdb database
    wget http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz && \\
        gzip -d VFDB_setA_nt.fas.gz && \\
        awk -v db=VFDB '/>/{ split(\$0,name," "); split(\$0,id," \\\\["); all=\$0; \$0=">" db "~~~" name[2] "~~~" name[1] "~~~[" id[2] "~~~" all }1' VFDB_setA_nt.fas | \\
        sed -e 's/~>/~/g' -e 's/ ~/~/g' -e 's/]~/~/g' -e 's/ >/ /' | \\
        awk -F "]" ' { if (\$0 ~ />/) { gsub(" ", "_", \$1); print \$1 "] " \$2 "]"} else { print \$0 }}' > sequences && \\
        makeblastdb -in sequences -title 'vfdb' -dbtype nucl -logfile /dev/null && \\
        rm VFDB_setA_nt.fas
    """
}
