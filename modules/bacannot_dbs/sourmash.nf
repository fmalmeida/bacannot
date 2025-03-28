process SOURMASH_DB {
    publishDir "${params.output}/sourmash_db", mode: 'copy', overwrite: "$params.force_update"
    label = [ 'process_ultralow', 'db_download' ]

    output:
    path "genbank-{21,31,51}.lca.json.gz"

    script:
    """
    # download sourmash database
    curl -L -o genbank-21.lca.json.gz https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k21.lca.json.gz
    curl -L -o genbank-31.lca.json.gz https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k31.lca.json.gz
    curl -L -o genbank-51.lca.json.gz https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k51.lca.json.gz
    """
}
