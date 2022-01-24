process PLASMIDFINDER {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename == "plasmidfinder") "plasmids/$filename"
    else null
  }
  tag "${prefix}"
  label 'main'

  input:
  tuple val(prefix), file(genome)

  output:
  tuple val(prefix), file("plasmidfinder") // Get everything
  tuple val(prefix), file("plasmidfinder/results_tab.tsv") // Get the main result

  script:
  """
  # Check thresholds
  [ "${params.plasmids_mincov}" > 1 ] && mincov=\$(( ${params.plasmids_mincov} / 100 )) || mincov=${params.plasmids_mincov}
  [ "${params.plasmids_minid}" > 1 ] && minid=\$(( ${params.plasmids_minid} / 100 )) || minid=${params.plasmids_minid}

  # Run plasmidfinder
  mkdir plasmidfinder ;
  plasmidfinder.py -i $genome -o plasmidfinder -l \$mincov -t \$minid -x

  # Remove tmp
  rm -rf plasmidfinder/tmp
  """
}
