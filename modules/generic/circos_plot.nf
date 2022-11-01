process CIRCOS {
    tag "$sample"

    label = [ 'misc', 'process_medium' ]

    input:
    tuple val(sample), path(assembly), path(inputs, stageAs: 'sample/*') // all inputs required by circos
    path(circos_conf, stageAs: 'conf/*')

    output:
    path("*")

    shell:
    '''
    # per contig
    while read contig ; do

        mkdir $contig      ;
        mkdir $contig/data ;
        cp -r conf $contig ;
        
        for file in $(ls sample/*) ; do

            name=$(basename $file) ;
            grep $contig $file > ${contig}/data/${name} || touch ${contig}/data/${name} ;
        
        done ;

        ( cd $contig && circos -conf conf/template.conf ) ;

    done < <( cat sample/skew.txt | cut -f 1 | sort -u )

    # normal one
    mkdir all             ;
    cp -r sample all/data ;
    cp -r conf all        ;
    ( cd all && circos -conf conf/template.conf ) ;
    '''
}