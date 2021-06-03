process victors {
  publishDir "${params.outdir}/${prefix}/virulence/victors", mode: 'copy'
  tag "Scanning virulence genes with Victors"
  label 'main'

  input:
  tuple val(prefix), file(genes)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_victors_blastp_onGenes.summary.txt")
  tuple val(prefix), file("${prefix}_victors_blastp_onGenes.txt")
  file('*.txt') // Grab summaries

  script:
  """
  # With predicted gene sequences

  run_blasts.py blastp --query $genes --db /work/dbs/victors/diamond.dmnd --minid ${params.blast_virulence_minid} \
  --mincov ${params.blast_virulence_mincov} --threads ${params.threads} --out ${prefix}_victors_blastp_onGenes.txt --2way | \
  sed -e 's/PRODUCT/VICTORS_ID/g' > ${prefix}_victors_blastp_onGenes.summary.txt ;
  """
}
