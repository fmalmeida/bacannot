process CIRCOS {
    publishDir "${params.output}/${prefix}/circos", mode: 'copy', saveAs: { filename ->
        if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
        else "$filename"
    }
    tag "$prefix"

    label = [ 'perl', 'process_low' ]

    input:
    tuple val(prefix), path(inputs)

    output:
    path("*")

    script:
    def genome = inputs[0]
    def gff    = inputs[1]

    """
    echo "${genome},${prefix},dgrey" > input.fofn

    plot_circos \\
        --fofn input.fofn \\
        --skip_links \\
        --bacannot \\
        --minlen ${params.circos_min_len} \\
        --outdir PLOT
    
    ( cd PLOT/conf && touch mges.txt )

    plot_circos --gff2labels "" "" ID black <( awk '\$7=="+"' ${gff} ) > PLOT/conf/forward_features.txt

    plot_circos --gff2labels "" "" ID dgreen <( awk '\$7=="-"' ${gff} ) > PLOT/conf/reverse_features.txt
    
    plot_circos --gff2labels "" "" ID dorange <( awk '\$3 ~ /rRNA/' ${gff} ) > PLOT/conf/rrna.txt
    
    plot_circos --gff2labels "" "" ID dpurple <( awk '\$3 ~ /tRNA/' ${gff} ) > PLOT/conf/trna.txt

    plot_circos --gff2labels "" "" "NDARO:Gene_Name" black <( awk '\$2 ~ /AMRFinderPlus/' ${gff} ) > PLOT/conf/bacannot_labels.txt

    plot_circos --gff2labels "" "" "VFDB:Product" black <( awk '\$2 ~ /VFDB/' ${gff} ) \\
        sed -e 's/[//g' -e 's/_(VF.* //g' > PLOT/conf/bacannot_labels.txt
    
    mv PLOT/* .
    rm -rf PLOT
    cd conf
    circos
    """
}