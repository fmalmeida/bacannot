process AMRFINDER2TSV {
    tag "$prefix"
    label = [ 'renv', 'process_low' ]

    input:
    tuple val(prefix), path(inputs)
    path(aro_tsv)

    output:
    tuple val(prefix), path("identified_amrfinderplus_genes.tsv"), emit: args, optional: true

    script:
    def annotation_summary = inputs[3]
    """
    amrfinder2tsv.R \
        --input $annotation_summary \
        --output identified_amrfinderplus_genes.tsv \
        --aro $aro_tsv \
        --sample $prefix
    """
}