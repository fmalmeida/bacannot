process SEQUENCESERVER {
    publishDir "${params.output}/${prefix}/SequenceServerDBs", mode: 'copy'
    tag "${prefix}"
    label = [ 'server', 'process_ultralow' ] 

    input:
    tuple val(prefix), file(genome), file(genes), file(proteins)

    output:
    tuple val(prefix), path("*")          , emit: all
    tuple val(prefix), path("${genome}")  , emit: genome
    tuple val(prefix), path("${genes}")   , emit: genes
    tuple val(prefix), path("${proteins}"), emit: proteins

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
