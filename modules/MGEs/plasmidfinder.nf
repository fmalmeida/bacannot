process PLASMIDFINDER {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename == "plasmidfinder") "plasmids/$filename"
    else null
  }
  tag "${prefix}"
  label = [ 'python', 'process_low' ]

  input:
  tuple val(prefix), file(genome)
  file(bacannot_db)

  output:
  tuple val(prefix), path("plasmidfinder")                , emit: all
  tuple val(prefix), path("plasmidfinder/results_tab.tsv"), emit: results

  script:
  """
  # Check thresholds
  [ "${params.plasmids_mincov}" > 1 ] && mincov=\$(( ${params.plasmids_mincov} / 100 )) || mincov=${params.plasmids_mincov}
  [ "${params.plasmids_minid}" > 1 ] && minid=\$(( ${params.plasmids_minid} / 100 )) || minid=${params.plasmids_minid}

  # Run plasmidfinder
  mkdir plasmidfinder ;
  plasmidfinder.py \\
      -i $genome \\
      -o plasmidfinder \\
      -l \$mincov \\
      -t \$minid \\
      -x \\
      --databasePath ${bacannot_db}/plasmidfinder_db

  # Remove tmp
  rm -rf plasmidfinder/tmp
  """
}
