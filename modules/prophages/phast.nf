process phast {
  publishDir "${params.outdir}/${prefix}/prophages/phast_db", mode: 'copy'
  tag "Annotating prophage genes with PHAST database"
  label 'main'

  input:
  tuple val(prefix), file(genes)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_phast_blastp_onGenes.summary.txt")
  tuple val(prefix), file("${prefix}_phast_blastp_onGenes.txt")
  file('*.txt') // Grab summaries

  script:
  """
  # With predicted gene sequences

  run_blasts.py blastp --query $genes --db /work/dbs/phast/diamond.dmnd --minid ${params.blast_MGEs_minid} \
  --mincov ${params.blast_MGEs_mincov} --threads ${params.threads} --out ${prefix}_phast_blastp_onGenes.txt --2way | \
  sed -e 's/PRODUCT/PHAST_ID/g' > ${prefix}_phast_blastp_onGenes.summary.txt ;
  """
}
