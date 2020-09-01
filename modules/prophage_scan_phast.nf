process phast {
  publishDir "${params.outdir}/${prefix}/prophages", mode: 'copy'
  tag "Scanning prophage genes with PHAST database"
  label 'main'

  input:
  tuple val(prefix), file(genes)
  tuple val(prefix), file(genome)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_phast_blastp_onGenes.summary.txt")
  tuple val(prefix), file("${prefix}_phast_blastx_onGenome.summary.txt")
  file('*.txt') // Grab summaries

  script:
  """
  # PHAST is a protein-only dabatase

  ## With predicted gene sequences
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastp --query $genes --db /work/dbs/phast/diamond.dmnd --minid ${params.diamond_MGEs_minid} \
  --mincov ${params.diamond_MGEs_mincov} --threads ${params.threads} --out ${prefix}_phast_blastp_onGenes.txt | \
  sed -e 's/PRODUCT/PHAST_ID/g' > ${prefix}_phast_blastp_onGenes.summary.txt ;

  ## With the whole genome
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastx --query $genome --db /work/dbs/phast/diamond.dmnd --minid ${params.diamond_MGEs_minid} \
  --mincov ${params.diamond_MGEs_mincov} --threads ${params.threads} --out ${prefix}_phast_blastx_onGenome.txt | \
  sed -e 's/PRODUCT/PHAST_ID/g' > ${prefix}_phast_blastx_onGenome.summary.txt ;
  """
}
