process argminer {
  publishDir "${params.outdir}/${prefix}/resistance/ARGMiner", mode: 'copy'
  tag "Scanning AMR genes with ARGMiner db"
  label 'main'

  input:
  tuple val(prefix), file(genes)
  tuple val(prefix), file(genome)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_argminer_blastp_onGenes.summary.txt")
  tuple val(prefix), file("${prefix}_argminer_blastx_onGenome.summary.txt")
  file('*.txt') // Grab summaries

  script:
  """
  # ARGMiner is a protein-only dabatase

  ## With predicted gene sequences
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastp --query $genes --db /work/dbs/argminer/diamond.dmnd --minid ${params.blast_resistance_minid} \
  --mincov ${params.blast_resistance_mincov} --threads ${params.threads} --out ${prefix}_argminer_blastp_onGenes.txt > ${prefix}_argminer_blastp_onGenes.summary.txt ;

  ## With the whole genome
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastx --query $genome --db /work/dbs/argminer/diamond.dmnd --minid ${params.blast_resistance_minid} \
  --mincov ${params.blast_resistance_mincov} --threads ${params.threads} --out ${prefix}_argminer_blastx_onGenome.txt > ${prefix}_argminer_blastx_onGenome.summary.txt ;
  """
}
