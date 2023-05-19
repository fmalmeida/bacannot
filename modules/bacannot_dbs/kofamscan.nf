process KOFAMSCAN_DB {
    publishDir "${params.output}/kofamscan_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_low' ]

    output:
    file("*")

    script:
    if (workflow.containerEngine != 'singularity') {
        chmod_cmd = 'chmod a+rw profiles.tar.gz ko_list'
        chown_cmd = 'chown -R root:\$(id -g) profiles'
        tar_cmd   = '--same-owner'
    } else {
        chmod_cmd = ''
        chown_cmd = ''
        tar_cmd   = ''
    }
    """
    # download kofamscan database
    wget --tries=10 ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
    wget --tries=10 ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
    gunzip ko_list.gz
    $chmod_cmd
    tar $tar_cmd -xvzf profiles.tar.gz
    $chown_cmd
    rm -rf profiles.tar.gz

    # for the sake of size and fastness
    # let's select only the KOs from prokaryotes
    cd profiles && \\
        for dirs in *.hmm ; do
            if ! grep -qxFe "\$dirs" prokaryote.hal ; then rm -rf \$dirs ; fi; 
        done
    """
}