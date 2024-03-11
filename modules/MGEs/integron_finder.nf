process INTEGRON_FINDER {
    publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
        if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
        else "integron_finder/$filename"
    }
    tag "${prefix}"
    label = [ 'process_medium' ]

    input:
    tuple val(prefix), file(genome)

    output:
    tuple val(prefix), path("*")                      , emit: all
    tuple val(prefix), path("${prefix}_integrons.gbk"), emit: gbk, optional: true
    path("integronfinder_version.txt")

    script:
    def args = task.ext.args ?: ''
    """
    # Get version
    integron_finder --version > integronfinder_version.txt ;

    # run tool
    integron_finder \\
        --local-max \\
        --func-annot \\
        --pdf \\
        --gbk \\
        --cpu $task.cpus \\
        $args \\
        $genome
    
    # move results
    mv Results_Integron_Finder_${prefix}/* . ;
    rm -rf Results_Integron_Finder_${prefix} ;
    
    # convert to gff if available
    for gbk in \$(ls *.gbk) ; do
        cat \$gbk >> ${prefix}_integrons.gbk ;
    done
    """
}
