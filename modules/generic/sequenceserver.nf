process sequenceserver {
    publishDir "${params.outdir}/${prefix}/SequenceServerDBs", mode: 'copy'
    tag "${prefix}"
    label 'main'

    input:
    tuple val(prefix), file(genome), file(genes), file(proteins)

    output:
    file("*")
    file("${genome}")
    file("${genes}")
    file("${proteins}")

    script:
    """
    # genome
    makeblastdb -in $genome -dbtype nucl -title "${prefix} genome" -parse_seqids

    # genes
    makeblastdb -in $genes -dbtype nucl -title "${prefix} genes" -parse_seqids

    # proteins
    makeblastdb -in $proteins -dbtype prot -title "${prefix} proteins" -parse_seqids
    """
}