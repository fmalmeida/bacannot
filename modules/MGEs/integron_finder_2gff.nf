process INTEGRON_FINDER_2GFF {
    publishDir "${params.output}/${prefix}/integron_finder", mode: 'copy'
    tag "${prefix}"
    label = [ 'misc', 'process_low' ]

    input:
    tuple val(prefix), file(gbk)

    output:
    tuple val(prefix), path("${prefix}_integrons.gff"), emit: gff

    script:
    def args = task.ext.args ?: ''
    """    
    # fix 0-based sequences
    sed -e 's/ 0\\.\\./ 1\\.\\./g' -e 's/complement(0\\.\\./complement(1\\.\\./g' $gbk > fixed.gbk
    
    # convert to gff if available
    conda run -n perl bp_genbank2gff3 fixed.gbk -o - | \\
        grep 'integron_id' | \\
        sed 's|ID=.*integron_id=|ID=|g' | \\
        sed 's/GenBank/Integron_Finder/g' >> ${prefix}_integrons.gff
    """
}
