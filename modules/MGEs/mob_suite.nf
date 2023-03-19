process MOBSUITE {
    publishDir "${params.output}", mode: 'copy', saveAs: { filename ->
        if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
        else "${prefix}/mob_suite/$filename"
    }
    tag "${prefix}"
    label = [ 'process_medium' ]

    conda "bioconda::mob_suite=3.1.4"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/mob_suite:3.1.4--pyhdfd78af_0' :
        'quay.io/biocontainers/mob_suite:3.1.4--pyhdfd78af_0' }"

    input:
    tuple val(prefix), file(genome)

    output:
    tuple val(prefix), path("${prefix}_mobtyper_results.txt"), emit: results
    path("mobtyper_version.txt")

    script:
    def args = task.ext.args ?: ''
    """
    # Get version
    mob_typer --version > mobtyper_version.txt ;

    # run tool
    mob_typer \\
        --multi \\
        --num_threads $task.cpus \\
        --sample_id $prefix \\
        --infile $genome \\
        $args \\
        --out_file ${prefix}_mobtyper_results.txt
    
    # convert to gff if available
    # for gbk in \$(ls *.gbk) ; do
    #     cat \$gbk >> ${prefix}_integrons.gbk ;
    # done
    """
}
