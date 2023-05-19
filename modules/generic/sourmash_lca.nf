process SOURMASH_LCA {
    publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
        if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
        else "sourmash/$filename"
    }
    tag "${prefix}"
    label = [ 'process_low' ]

    input:
    path bacannot_dbs
    tuple val(prefix), path(genome)
    val scale
    val kmer

    output:
    path "*"
    path "sourmash_version.txt"

    when: !params.skip_sourmash

    script:
    def lca_db = "${bacannot_dbs}/sourmash_db/genbank-k${kmer}.lca.json.gz"
    """
    # get version file
    sourmash --version > sourmash_version.txt

    # sketch input genome
    sourmash \\
        sketch dna \\
        -p scaled=${scale},k=${kmer} \\
        --name-from-first \\
        ${genome}

    # classify
    sourmash \\
        lca classify \\
        --db ${lca_db} \\
        --query ${genome}.sig
    
    # summarize
    sourmash \\
        lca summarize \\
        --db ${lca_db} \\
        --query some-genome.fa.gz.sig > ${genome}_sourmash.summary.txt
    """
}
