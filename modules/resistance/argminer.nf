process ARGMINER {
  publishDir "${params.output}/${prefix}/resistance/ARGMiner", mode: 'copy'
  tag "${prefix}"
  label = [ 'misc', 'process_low' ]

  input:
  tuple val(prefix), file(genes)
  file(bacannot_db)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_argminer_blastp_onGenes.summary.txt")
  file('*.txt') // Grab summaries

  script:
  """
  # run blast with argminer db
  run_blasts.py \\
      blastp \\
      --query $genes \\
      --db ${bacannot_db}/argminer_db/diamond.dmnd \\
      --minid ${params.blast_resistance_minid} \\
      --mincov ${params.blast_resistance_mincov} \\
      --threads $task.cpus \\
      --out ${prefix}_argminer_blastp_onGenes.txt \\
      --2way > ${prefix}_argminer_blastp_onGenes.summary.txt ;
  """
}
