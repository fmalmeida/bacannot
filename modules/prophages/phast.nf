process PHAST {
  publishDir "${params.output}/${prefix}/prophages/phast_db", mode: 'copy'
  tag "${prefix}"
  label = [ 'misc', 'process_low' ]

  input:
  tuple val(prefix), file(genes)
  file(bacannot_db)

  output:
  tuple val(prefix), path("${prefix}_phast_blastp_onGenes.summary.txt"), emit: summary
  tuple val(prefix), path("${prefix}_phast_blastp_onGenes.txt")        , emit: results
  tuple val(prefix), path('*.txt')                                     , emit: all

  script:
  """
  # PHAST has protein database
  run_blasts.py \\
      blastp \\
      --query $genes \\
      --db ${bacannot_db}/phast_db/diamond.dmnd \\
      --minid ${params.blast_MGEs_minid} \\
      --mincov ${params.blast_MGEs_mincov} \\
      --threads $task.cpus \\
      --out ${prefix}_phast_blastp_onGenes.txt --2way | \\
  sed -e 's/PRODUCT/PHAST_ID/g' -e 's/;//g' > ${prefix}_phast_blastp_onGenes.summary.txt ;
  """
}
