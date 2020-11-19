process iceberg {
  publishDir "${params.outdir}/${prefix}/ICEs", mode: 'copy'
  tag "Scanning for ICE genes with ICEberg database"
  label 'main'

  input:
  tuple val(prefix), file(genes_aa)
  tuple val(prefix), file(genome)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_iceberg_blastp_onGenes.summary.txt")
  tuple val(prefix), file("${prefix}_iceberg_blastp_onGenes.txt")
  tuple val(prefix), file("${prefix}_iceberg_blastn_onGenome.summary.txt")
  file('*.txt') // Grab all

  script:
  """
  # ICEberg is a protein and nucleotide dabatase
  # In protein are the genes found inside ICEs
  # In nucleotide are the full-length ICEs

  ## Checking ICE genes

  ## With predicted gene sequences
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastp --query $genes_aa --db /work/dbs/iceberg/diamond.dmnd --minid ${params.blast_MGEs_minid} \
  --mincov ${params.blast_MGEs_mincov} --threads ${params.threads} --out ${prefix}_iceberg_blastp_onGenes.txt --2way | \
  sed -e 's/GENE/ICEBERG_ID/g' > ${prefix}_iceberg_blastp_onGenes.summary.txt ;

  ## Checking for full-length ICEs
  ### The blast db was throwing errors
  makeblastdb -dbtype nucl -in /work/dbs/iceberg/sequences -out sequences ;
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastn --query $genome --db sequences --minid 0 --mincov 0 --threads ${params.threads} \
  --out ${prefix}_iceberg_blastn_onGenome.txt | sed -e 's/GENE/ICEBERG_ID/g' > ${prefix}_iceberg_blastn_onGenome.summary.txt ;
  """
}
