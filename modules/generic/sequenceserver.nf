process SEQUENCESERVER {
    publishDir "${params.output}/${prefix}/SequenceServerDBs", mode: 'copy'
    tag "${prefix}"
    label = [ 'server', 'process_ultralow' ]
    

    input:
    tuple val(prefix), file(genome), file(genes), file(proteins)

    output:
    path("*")          , emit: all
    path("${genome}")  , emit: genome
    path("${genes}")   , emit: genes
    path("${proteins}"), emit: proteins

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
