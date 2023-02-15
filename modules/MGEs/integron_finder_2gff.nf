process INTEGRON_FINDER_2GFF {
    publishDir "${params.output}", mode: 'copy'
    tag "${prefix}"
    label = [ 'misc', 'process_low' ]

    input:
    tuple val(prefix), file(gbk)

    output:
    tuple val(prefix), path("${prefix}_integrons.gff"), emit: gff, optional: true

    script:
    def args = task.ext.args ?: ''
    """    
    # convert to gff if available
    touch ${prefix}_integrons.gff ;
    for gbk in \$(ls *.gbk) ; do
        conda run -n perl bp_genbank2gff3 \$gbk -o - | \
            grep 'integron_id' | \
            sed 's|ID=.*integron_id=|ID=|g' >> ${prefix}_integrons.gff
    done
    """
}
