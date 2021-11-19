process AMRFINDER_DB {
    publishDir "${params.output}/amrfinder_db", mode: 'copy', overwrite: "$params.force_update"
    label 'db_download'
   
    output:
    file("*")

    script:
    """
    # download amrfinderplus database
    wget -rkN -e robots=off --no-parent --reject "index.html*" https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/

    # save only database files
    mv ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/* .
    rm -r ftp.ncbi.nlm.nih.gov
    """
}
