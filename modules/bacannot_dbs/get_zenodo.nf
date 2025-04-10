process GET_ZENODO_DB {
    publishDir "${params.output}", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'db_download', 'process_low' ]

    tag "Downloading pre-built databases"

    output:
    file("*")

    script:
    """
    # download database from zenodo
    zenodo_get 12674473

    # organize data
    tar zxvf *.tar.gz && rm *.tar.gz
    rm -rf \$( find . -name 'pipeline_info' )
    """
}