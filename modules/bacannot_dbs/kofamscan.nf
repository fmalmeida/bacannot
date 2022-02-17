process KOFAMSCAN_DB {
    publishDir "${params.output}/kofamscan_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_low' ]
   
    output:
    file("*")

    script:
    """
    # download kofamscan database
    wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz && \\
        wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz && \\
        gunzip ko_list.gz && \\
        tar xvzf profiles.tar.gz && \\
        rm -rf profiles.tar.gz

    # for the sake of size and fastness
    # let's select only the KOs from prokaryotes
    cd profiles && \\
        for dirs in *.hmm ; do
            if ! grep -qxFe "\$dirs" prokaryote.hal ; then rm -rf \$dirs ; fi; 
        done
    """
}