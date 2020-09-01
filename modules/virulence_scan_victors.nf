process victors {
  publishDir "${params.outdir}/${prefix}/virulence/Victors", mode: 'copy'
  tag "Scanning virulence genes with Victors"
  label 'main'

  input:
  tuple val(prefix), file(genes)
  tuple val(prefix), file(genome)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_victors_blastp_onGenes.summary.txt")
  tuple val(prefix), file("${prefix}_victors_blastx_onGenome.summary.txt")
  file('*.txt') // Grab summaries

  script:
  """
  # Victors is a protein-only dabatase

  ## With predicted gene sequences
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastp --query $genes --db /work/dbs/victors/diamond.dmnd --minid ${params.diamond_virulence_minid} \
  --mincov ${params.diamond_virulence_mincov} --threads ${params.threads} --out ${prefix}_victors_blastp_onGenes.txt | \
  sed -e 's/PRODUCT/VICTORS_ID/g' > ${prefix}_victors_blastp_onGenes.summary.txt ;

  ## With the whole genome
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastx --query $genome --db /work/dbs/victors/diamond.dmnd --minid ${params.diamond_virulence_minid} \
  --mincov ${params.diamond_virulence_mincov} --threads ${params.threads} --out ${prefix}_victors_blastx_onGenome.txt | \
  sed -e 's/PRODUCT/VICTORS_ID/g' > ${prefix}_victors_blastx_onGenome.summary.txt ;
  """
}
