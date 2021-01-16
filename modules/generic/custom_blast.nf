process custom_blast {
  publishDir "${params.outdir}/${prefix}/custom_annotations/${customDB.baseName}", mode: 'copy'
  tag "Performing annotation with User's custom db"
  label 'main'

  input:
  tuple val(prefix), file(genome)
  each file(customDB)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_${customDB.baseName}_blastn.summary.txt")
  file('*.txt') // Grab all

  script:
  """
  # Step 1 - Create blast db
  makeblastdb -in $customDB -dbtype nucl -out customDB -logfile /dev/null ;

  # Step 2 - Execute blastn
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastn --query $genome --db customDB --minid ${params.blast_custom_minid} \
  --mincov ${params.blast_custom_mincov} --threads ${params.threads} --out ${prefix}_${customDB.baseName}_blastn.txt > ${prefix}_${customDB.baseName}_blastn.summary.txt ;
  """
}
