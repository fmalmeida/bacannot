process SOURMASH_ALL {
    publishDir "${params.output}/sourmash_all", mode: 'copy', saveAs: { filename ->
        if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
        else "$filename"
    }
    label = [ 'process_medium' ]

    input:
    path( "genomes/*" )
    val scale
    val kmer

    output:
    path "*"                       , emit: all
    path "sourmash_version.txt"    , emit: versions
    path "sourmash_cmp.matrix.png" , emit: plot

    when: !params.skip_sourmash

    script:
    """
    # get version file
    sourmash --version > sourmash_version.txt

    # sketch input genomes
    mkdir signatures ;
    ( 
        cd genomes && \\
        for genome in * ; do
            export name=\$( echo \${genome} | cut -f 2 -d '/' | cut -f 1,2 -d '_' ) ; 
            sourmash \\
                sketch dna \\
                -p scaled=${scale},k=${kmer} \\
                \${genome} \\
                -o ../signatures/\${name}.sig ;
        done ;
    )

    # compare
    sourmash \\
        compare \\
        signatures/* \\
        -p $task.cpus \\
        -o sourmash_cmp
    
    # plot
    sourmash plot --labels sourmash_cmp
    sourmash \\
        plot \\
        --pdf \\
        --csv sourmash_plot.csv \\
        --labels sourmash_cmp
    """
}
