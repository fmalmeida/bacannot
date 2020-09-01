process vfdb {
  publishDir "${params.outdir}/${prefix}/virulence/VFDB", mode: 'copy'
  tag "Scanning virulence genes with VFDB"
  label 'main'

  input:
  tuple val(prefix), file(genes)
  tuple val(prefix), file(genome)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_vfdb_blastn_onGenes.summary.txt")
  tuple val(prefix), file("${prefix}_vfdb_blastn_onGenome.summary.txt")
  file('*.txt') // Grab summaries

  script:
  """
  # VFDB is a nucleotide-only dabatase

  ## With predicted gene sequences
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastn --query $genes --db /work/dbs/vfdb/sequences --minid ${params.diamond_virulence_minid} \
  --mincov ${params.diamond_virulence_mincov} --threads ${params.threads} --out ${prefix}_vfdb_blastn_onGenes.txt | \
  sed -e 's/ACCESSION/VFDB_ID/g' > ${prefix}_vfdb_blastn_onGenes.summary.txt ;

  ## With the whole genome
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastn --query $genome --db /work/dbs/vfdb/sequences --minid ${params.diamond_virulence_minid} \
  --mincov ${params.diamond_virulence_mincov} --threads ${params.threads} --out ${prefix}_vfdb_blastn_onGenes.txt | \
  sed -e 's/ACCESSION/VFDB_ID/g' > ${prefix}_vfdb_blastn_onGenome.summary.txt ;
  """

}
