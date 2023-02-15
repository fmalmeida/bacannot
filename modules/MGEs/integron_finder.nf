process INTEGRON_FINDER {
    publishDir "${params.output}", mode: 'copy', saveAs: { filename ->
        if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
        else "${prefix}/integron_finder/$filename"
    }
    tag "${prefix}"
    label = [ 'process_medium' ]

    conda "bioconda::integron_finder=2.0.1"
        container "${ workflow.containerEngine == 'singularity' ?
            'https://depot.galaxyproject.org/singularity/integron_finder:2.0.1--pyhdfd78af_0' :
            'quay.io/biocontainers/integron_finder:2.0.1--pyhdfd78af_0' }"

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
