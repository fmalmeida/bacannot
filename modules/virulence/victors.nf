process VICTORS {
  publishDir "${params.output}/${prefix}/virulence/victors", mode: 'copy'
  tag "${prefix}"
  label = [ 'misc', 'process_low' ]

  input:
  tuple val(prefix), file(genes)
  file(bacannot_db)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), path("${prefix}_victors_blastp_onGenes.summary.txt")
  tuple val(prefix), path("${prefix}_victors_blastp_onGenes.txt")
  path('*.txt')

  script:
  """
  # Victors has protein database
  run_blasts.py \\
      blastp \\
      --query $genes \\
      --db ${bacannot_db}/victors_db/diamond.dmnd \\
      --minid ${params.blast_virulence_minid} \\
      --mincov ${params.blast_virulence_mincov} \\
      --threads ${params.threads} \\
      --out ${prefix}_victors_blastp_onGenes.txt \\
      --2way | \\
  sed -e 's/PRODUCT/VICTORS_ID/g' > ${prefix}_victors_blastp_onGenes.summary.txt ;
  """
}
