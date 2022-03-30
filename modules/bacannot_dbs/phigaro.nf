process PHIGARO_DB {
    publishDir "${params.output}/phigaro_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_medium' ]
   
    output:
    file("*")

    script:
    """   
    # download phigaro database
    wget --tries=10 http://download.ripcm.com/phigaro/allpvoghmms
    wget --tries=10 http://download.ripcm.com/phigaro/allpvoghmms.h3f
    wget --tries=10 http://download.ripcm.com/phigaro/allpvoghmms.h3i
    wget --tries=10 http://download.ripcm.com/phigaro/allpvoghmms.h3m
    wget --tries=10 http://download.ripcm.com/phigaro/allpvoghmms.h3p
    """
}
