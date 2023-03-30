process MLST {
    publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
        if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
        else "MLST/$filename"
    }
    tag "${prefix}"
    label = [ 'process_ultralow' ]

    input:
    tuple val(prefix), file(genome)
    file(bacannot_db)

    output:
    tuple val(prefix), path("*")                            , emit: all
    tuple val(prefix), path("${prefix}_mlst_analysis.txt")  , emit: results optional true
    tuple val(prefix), path("${prefix}_novel_alleles.fasta"), emit: alleles optional true
    path('mlst_version.txt')                                , emit: version

    script:
    """
    # update tool database
    mlst-make_blast_db.sh ${bacannot_db}/mlst_db

    # Save mlst tool version
    mlst --version > mlst_version.txt ;

    # run mlst
    mlst \\
        --quiet \\
        --novel ${prefix}_novel_alleles.fasta \\
        $genome > ${prefix}_mlst_analysis.txt
    """
}
