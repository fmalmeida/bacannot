process VFDB2TSV {
    tag "$prefix"
    label = [ 'renv', 'process_low' ]

    input:
    tuple val(prefix), path(inputs)

    output:
    tuple val(prefix), path("identified_vfdb_genes.tsv"), emit: genes, optional: true

    script:
    def annotation_summary = inputs[3]
    """
    vfdb2tsv.R \
        --input $annotation_summary \
        --output identified_vfdb_genes.tsv \
        --sample $prefix
    """
}