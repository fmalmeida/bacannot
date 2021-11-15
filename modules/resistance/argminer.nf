process argminer {
  publishDir "${params.outdir}/${prefix}/resistance/ARGMiner", mode: 'copy'
  tag "${prefix}"
  label 'main'

  input:
  tuple val(prefix), file(genes)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_argminer_blastp_onGenes.summary.txt")
  file('*.txt') // Grab summaries

  script:
  """
  # With predicted gene sequences

  run_blasts.py blastp --query $genes --db /work/dbs/ARGMiner/diamond.dmnd --minid ${params.blast_resistance_minid} \
  --mincov ${params.blast_resistance_mincov} --threads ${params.threads} --out ${prefix}_argminer_blastp_onGenes.txt --2way > ${prefix}_argminer_blastp_onGenes.summary.txt ;
  """
}
