process vfdb {
  publishDir "${params.outdir}/${prefix}/virulence/VFDB", mode: 'copy'
  tag "Scanning virulence genes with VFDB"
  label 'main'

  input:
  tuple val(prefix), file(genes)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_vfdb_blastn_onGenes.summary.txt")
  tuple val(prefix), file("${prefix}_vfdb_blastn_onGenes.txt")
  file('*.txt') // Grab summaries

  script:
  """
  # With predicted gene sequences
  
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastn --query $genes --db /work/dbs/vfdb/sequences --minid ${params.blast_virulence_minid} \
  --mincov ${params.blast_virulence_mincov} --threads ${params.threads} --out ${prefix}_vfdb_blastn_onGenes.txt --2way | \
  sed -e 's/ACCESSION/VFDB_ID/g' > ${prefix}_vfdb_blastn_onGenes.summary.txt ;
  """

}
