process MERGE_SUMMARIES {
    publishDir "${params.output}", mode: 'copy'
    label = [ 'misc', 'process_low' ]

    input:
    path(summaries)

    output:
    path("bacannot_summary.json"), emit: summary

    script:
    '''
    jq -s 'reduce .[] as $item ({}; . * $item)' *.json > bacannot_summary.json
    '''
}
