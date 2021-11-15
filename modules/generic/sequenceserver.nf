process sequenceserver {
    publishDir "${params.output}/${prefix}/SequenceServerDBs", mode: 'copy'
    tag "Generate SequenceServer DBs"
    label 'main'

    input:
    tuple val(prefix), file(genome), file(genes), file(proteins)

    output:
    file("*")

    script:
    """
    # genome
    makeblastdb -in $genome -dbtype nucl -title "${prefix} genome" -out ${prefix}_genome

    # genes
    makeblastdb -in $genes -dbtype nucl -title "${prefix} genes" -out ${prefix}_genes

    # proteins
    makeblastdb -in $proteins -dbtype prot -title "${prefix} proteins" -out ${prefix}_proteins
    """
}