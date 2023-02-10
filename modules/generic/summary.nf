process SUMMARY {
    publishDir "${params.output}/${prefix}", mode: 'copy'
    tag "${prefix}"
    label = [ 'misc', 'process_low' ]
    

    input:
    tuple val(prefix), 
    file(annotation), file(stageAs: "results/${prefix}/MLST/*"), 
    file(stageAs: "results/${prefix}/rRNA/*"), file(stageAs: "results/${prefix}/*"), 
    file(stageAs: "results/${prefix}/plasmids/*"), file(stageAs: "results/${prefix}/plasmids/*"), 
    file(stageAs: "results/${prefix}/genomic_islands/*"), file(stageAs: "results/${prefix}/virulence/vfdb/*"),
    file(stageAs: "results/${prefix}/virulence/victors/*"), file(stageAs: "results/${prefix}/prophages/phast_db/*"),
    file(stageAs: "results/${prefix}/prophages/phigaro/*"), file(stageAs: "results/${prefix}/prophages/*"),
    file(stageAs: "results/${prefix}/ICEs/*"), file(stageAs: "results/${prefix}/resistance/AMRFinderPlus/*"),
    file(stageAs: "results/${prefix}/resistance/RGI/*"), file(stageAs: "results/${prefix}/resistance/ARGMiner/*"),
    file(stageAs: "results/${prefix}/resistance/*"), file(stageAs: "results/${prefix}/methylations/*"),
    file(stageAs: "results/${prefix}/refseq_masher/*"), file(stageAs: "results/${prefix}/*"),
    file(stageAs: "results/${prefix}/*"), file(stageAs: "results/${prefix}/gffs/*")

    output:
    tuple val(prefix), path("${prefix}_summary.json"), emit: summaries

    script:
    """
    mkdir -p results/${prefix}/annotation
    ln -rs annotation/* results/${prefix}/annotation
    source activate falmeida-py
    falmeida-py bacannot2json -i results -o ${prefix}_summary.json
    """
}
