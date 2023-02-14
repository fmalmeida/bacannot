process INTEGRON_FINDER_2GFF {
    publishDir "${params.output}", mode: 'copy'
    tag "${prefix}"
    label = [ 'misc', 'process_low' ]

    input:
    tuple val(prefix), file(gbk)

    output:
    tuple val(prefix), path("${meta.id}_integrons.gff"), emit: gff, optional: true

    script:
    def args = task.ext.args ?: ''
    """
    # convert to gff if available
    touch ${meta.id}_integrons.gff ;
    for gbk in \$(ls *.gbk) ; do
        bp_genbank2gff3 $integron_finder -o - | \
            grep 'integron_id' | \
            sed 's|ID=.*integron_id=|ID=|g' >> ${meta.id}_integrons.gff
    done
    """
}
