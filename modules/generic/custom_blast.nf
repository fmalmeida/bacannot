process custom_blast {
  publishDir "${params.outdir}/${prefix}/custom_annotations/${customDB.baseName}", mode: 'copy'
  tag "Performing annotation with User's custom db"
  label 'main'

  input:
  tuple val(prefix), file(gff), file(genome)
  each file(customDB)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), val("${customDB.baseName}"), file("${prefix}_${customDB.baseName}_blastn.txt"), file("${prefix}_${customDB.baseName}_blastn.gff")
  file('*.txt') // Grab all

  script:
  """
  # Step 1 - Create blast db
  makeblastdb -in $customDB -dbtype nucl -out customDB ;

  # Step 2 - Execute blastn
  /miniconda/bin/python3 /usr/local/bin/run_blasts.py blastn --query $genome --db customDB --minid ${params.blast_custom_minid} \
  --mincov ${params.blast_custom_mincov} --threads ${params.threads} --out ${prefix}_${customDB.baseName}_blastn.txt > ${prefix}_${customDB.baseName}_blastn.summary.txt ;

  # Step 3 - Get BED from blastn
  awk '{print \$1 "\t" \$2 "\t" \$3}' ${prefix}_${customDB.baseName}_blastn.txt | tail -n +2 > ${prefix}_${customDB.baseName}_blastn.bed ;

  # Step 4 - Find intersection with annotation
  bedtools intersect -wa -a $gff -b ${prefix}_${customDB.baseName}_blastn.bed > ${prefix}_${customDB.baseName}_blastn.gff ;
  """
}
