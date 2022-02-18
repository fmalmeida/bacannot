process PHAST {
  publishDir "${params.output}/${prefix}/prophages/phast_db", mode: 'copy'
  tag "${prefix}"
  label = [ 'misc', 'process_low' ]

  input:
  tuple val(prefix), file(genes)
  file(bacannot_db)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), path("${prefix}_phast_blastp_onGenes.summary.txt")
  tuple val(prefix), path("${prefix}_phast_blastp_onGenes.txt")
  path('*.txt')

  script:
  """
  # PHAST has protein database
  run_blasts.py \\
      blastp \\
      --query $genes \\
      --db ${bacannot_db}/phast_db/diamond.dmnd \\
      --minid ${params.blast_MGEs_minid} \\
      --mincov ${params.blast_MGEs_mincov} \\
      --threads ${params.threads} \\
      --out ${prefix}_phast_blastp_onGenes.txt --2way | \\
  sed -e 's/PRODUCT/PHAST_ID/g' > ${prefix}_phast_blastp_onGenes.summary.txt ;
  """
}
