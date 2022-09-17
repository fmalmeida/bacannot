process VICTORS {
  publishDir "${params.output}/${prefix}/virulence/victors", mode: 'copy'
  tag "${prefix}"
  label = [ 'misc', 'process_low' ]

  input:
  tuple val(prefix), file(genes)
  file(bacannot_db)

  output:
  tuple val(prefix), path("${prefix}_victors_blastp_onGenes.summary.txt"), emit: summary
  tuple val(prefix), path("${prefix}_victors_blastp_onGenes.txt")        , emit: results
  path('*.txt')                                                          , emit: all

  script:
  """
  # Victors has protein database
  run_blasts.py \\
      blastp \\
      --query $genes \\
      --db ${bacannot_db}/victors_db/diamond.dmnd \\
      --minid ${params.blast_virulence_minid} \\
      --mincov ${params.blast_virulence_mincov} \\
      --threads $task.cpus \\
      --out ${prefix}_victors_blastp_onGenes.txt \\
      --2way | \\
  sed -e 's/PRODUCT/VICTORS_ID/g' > ${prefix}_victors_blastp_onGenes.summary.txt ;
  """
}
