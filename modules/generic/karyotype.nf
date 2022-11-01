process MAKE_KARYOTYPE {
    tag "$prefix"

    input:
    tuple val(prefix), path(inputs)

    output:
    tuple val(prefix), path("karyotype.txt"), emit: karyotype

    shell:
    '''
    while read -r FASTA FASTA_PREFIX FASTA_COLOR ; do

        name="$(basename $FASTA)" ;

        [ -s !{inputs[0]} ] && \
        bioawk \
            -c fastx \
            -v p=$FASTA_PREFIX \
            -v color=$FASTA_COLOR \
            '{ printf "chr - " substr($name,1) " " substr($name,1) " " "0" " " length($seq) " " color"\\n" }' \
            $FASTA >> karyotype.txt \
        || true ;
    
    done < <(echo "!{inputs[0]}\t!{prefix}\t!{params.karyotype_color}")
    '''
}