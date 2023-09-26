process GC_SKEW {
    tag "$prefix"
    label = [ 'misc', 'process_low' ]

    input:
    tuple val(prefix), path(inputs)

    output:
    tuple val(prefix), path('skew.txt'), emit: skew

    shell:
    '''
    if [ -s !{inputs[0]} ]; then

    # split per contig
    cat !{inputs[0]} | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".splitted.fasta")}
        print $0 >> filename
        close(filename)
    }'

    # iterate on each and get best size of window
    for contig in $(ls *.splitted.fasta) ; do
        seq_len=$(grep -v '>' $contig | awk '{ c+=length($0);} END { print c; }' | tr -d '' | tr -d '\n') ;

        if [ $seq_len -le 100000 ] ; then
            window=1000 ;
        elif [ $seq_len -le 200000 ] ; then
            window=5000 ;
        else
            window=20000 ;
        fi

        # GCcalc.py -f ${RESULTS}/all_vs_all_blast/concatenated_genomes.fasta -w $GCWINDOW -s $GCSTEP | \
    cut -f 1,2,3,5 | awk '{ if ($4 > 0) print $0 "\t" "color=dblue"; else print $0 "\t" "color=red"}' > ${RESULTS}/conf/GC_skew.txt

        GCcalc.py \
            -f $contig \
            -w $window \
            -s $window | \
            cut -f 1,2,3,5 | \
            awk '{ if ($4 > 0) print $0 "\t" "color=black"; else print $0 "\t" "color=color6" }' >> skew.txt ;
    
    done

    fi
    '''
}