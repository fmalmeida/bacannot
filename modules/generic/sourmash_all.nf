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
    path "*"
    path "sourmash_version.txt"

    when: !params.skip_sourmash

    script:
    """
    # get version file
    sourmash --version > sourmash_version.txt

    # sketch input genomes
    mkdir signatures ;
    for genome in genomes/* ; do
        export name=\$( echo \${genome} | cut -f 2 -d '/' | cut -f 1,2 -d '_' ) ; 
        sourmash \\
            sketch dna \\
            -p scaled=${scale},k=${kmer} \\
            \${genome} \\
            -o signatures/\${name}.sig ;
    done

    # compare
    sourmash \\
        compare \\
        signatures/* \\
        -p $task.cpus \\
        -o sourmash_cmp
    
    # plot
    sourmash \\
        plot \\
        --pdf \\
        --labels sourmash_cmp
    """
}
