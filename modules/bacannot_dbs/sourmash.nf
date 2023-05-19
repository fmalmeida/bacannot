process SOURMASH_DB {
    publishDir "${params.output}/sourmash_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'process_ultralow', 'db_download' ]

    output:
    path "genbank-{21,31,51}.lca.json.gz"

    script:
    """
    # download sourmash database
    curl -L -o genbank-21.lca.json.gz https://osf.io/gk2za/download
    curl -L -o genbank-31.lca.json.gz https://osf.io/ypsjq/download
    curl -L -o genbank-51.lca.json.gz https://osf.io/297dp/download
    """
}
