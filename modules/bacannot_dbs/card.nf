process CARD_DB {
    publishDir "${params.output}/card_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_ultralow' ]

    output:
    file("*")

    script:
    """
    # download CARD database
    wget --tries=10 https://card.mcmaster.ca/latest/data
    tar -xvf data ./card.json
    rm data
    """
}
