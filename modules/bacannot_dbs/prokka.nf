process PROKKA_DB {
    publishDir "${params.output}/prokka_db", mode: 'copy', overwrite: "$params.force_update"
    label 'db_download'
   
    output:
    file("*")

    script:
    """
    # download prokka additional database
    wget https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/TIGRFAMs_15.0_HMM.LIB.gz && \
	    gzip -d TIGRFAMs_15.0_HMM.LIB.gz
    """
}
