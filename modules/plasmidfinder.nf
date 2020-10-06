process plasmidfinder {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  tag "Detecting plasmids with plasmidfinder"
  label 'main'

  input:
  tuple val(prefix), file(genome)

  output:
  tuple val(prefix), file("plasmidfinder_results") // Get everything

  script:
  if ("${params.plasmids_mincov}" > 1) {
    params.plasmids_mincov = params.plasmids_mincov / 100
  }

  if (params.plasmids_minid > 1) {
    params.plasmids_minid = params.plasmids_mincov / 100
  }
  """
  # Check thresholds
  [ "${params.plasmids_mincov}" > 1 ] && mincov=\$(( ${params.plasmids_mincov} / 100 )) || mincov=${params.plasmids_mincov}
  [ "${params.plasmids_minid}" > 1 ] && minid=\$(( ${params.plasmids_minid} / 100 )) || minid=${params.plasmids_minid}

  # Activate conda environment
  source activate PLASMIDFINDER ;

  # Run plasmidfinder
  mkdir plasmidfinder_results ;
  plasmidfinder.py -i $genome -o plasmidfinder_results -l \$mincov -t \$minid -x
  """
}
