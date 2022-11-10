process PLOT {
    publishDir "${params.output}/${prefix}/CIRCOS", mode: 'copy', saveAs: { filename ->
        if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
        else "$filename"
    }
    tag "$prefix"

    label = [ 'perl', 'process_low' ]

    input:
    tuple val(prefix), path(inputs, stageAs: 'sample/*') // all inputs required by circos
    path(circos_conf, stageAs: 'conf/*')

    output:
    path("per_contig")
    path("whole_genome")

    shell:
    '''
    # per contig
    while read contig ; do

        mkdir -p per_contig/${contig}   ;
        mkdir per_contig/${contig}/data ;
        cp -r conf per_contig/${contig} ;
        
        for file in $(ls sample/*) ; do

            name=$(basename $file) ;
            grep $contig $file > per_contig/${contig}/data/${name} || touch per_contig/${contig}/data/${name} ;
        
        done ;

        ( cd per_contig/${contig} && circos -conf conf/template.conf ) ;

    done < <( cat sample/skew.txt | cut -f 1 | sort -u )

    # normal one
    mkdir whole_genome             ;
    cp -r sample whole_genome/data ;
    cp -r conf whole_genome        ;
    ( cd whole_genome && circos -conf conf/template.conf ) ;
    '''
}