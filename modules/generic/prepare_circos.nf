process PREPARE_CIRCOS {
    tag "$prefix"

    label = [ 'misc', 'process_low' ]

    input:
    tuple val(prefix), path(inputs, stageAs: '*') // all inputs required by circos

    output:
    tuple val(prefix), path("*.txt"), emit: data

    shell:
    '''
    # prepare amrfinderplus genes
    touch amrfinder.txt amrfinder_text.txt
    if [ -s identified_amrfinderplus_genes.tsv ]; then
        tail -n+2 \\
            identified_amrfinderplus_genes.tsv | \\
            cut -f 2,3,4 > amrfinder.txt ;
        tail -n+2 \\
            identified_amrfinderplus_genes.tsv | \\
            cut -f 2,3,4,5 > amrfinder_text.txt ;
    fi

    # prepare vfdb genes
    touch vfdb.txt vfdb_text.txt
    if [ -s identified_vfdb_genes.tsv ]; then
        tail -n+2 \\
            identified_vfdb_genes.tsv | \\
            cut -f 2,3,4 > vfdb.txt ;
        tail -n+2 \\
            identified_vfdb_genes.tsv | \\
            cut -f 2,3,4,5 > vfdb_text.txt ;
    fi

    # prepare plasmidfinder loci
    touch plasmidfinder.txt plasmidfinder_text.txt
    find \\
        -L . \\
        -name 'results_tab.tsv' | \\
        grep plasmidfinder | \\
        grep !{prefix} > tmpFile || true
    if [ -s tmpFile ] ; then

        plasmidfinder_results=$(cat tmpFile | tr -d ' ' | tr -d '\\n') ;
        if [ -s $plasmidfinder_results ]; then
            
            awk '
                BEGIN { OFS="\\t";}
                !/Database/ { print $7,$8,$2 }
            ' $plasmidfinder_results | \\
            tr '.' '_' > plasmidfinder_text.txt
            cut -f 1,2 plasmidfinder_text.txt > plasmidfinder.txt
            sed -i "s/__/\\t/g" plasmidfinder*txt

        fi

    fi

    # collapse labels
    cat *_text.txt > all_labels.txt
    '''
}