process CUSTOM_DATABASE {
  publishDir "${params.output}/${prefix}/custom_annotations/${customDB.baseName}", mode: 'copy'
  tag "${prefix}"
  label 'main'

  input:
  tuple val(prefix), file(gff), file(genome), file(proteins)
  each file(customDB)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), val("${customDB.baseName}"), file("${prefix}_${customDB.baseName}*.summary.txt"), file("${prefix}_${customDB.baseName}*.gff")
  file('*.txt') // Grab all

  script:
  """
  # Step 1 - Check if input is nucl or protein
  if seqkit stats input.faa | grep -iq "Protein" ; then
      blast_db=prot ; blast_cmd="blastp" ; blast_subj=${proteins} ;
  else
      blast_db=nucl ; blast_cmd="blastn" ; blast_subj=${genome} ;
  fi

  # Step 2 - Create blast db
  makeblastdb -in $customDB -dbtype \$blast_db -out customDB ;

  # Step 3 - Execute blast
  run_blasts.py \\
      \$blast_cmd \\
      --query \$blast_subj \\
      --db customDB \\
      --minid ${params.blast_custom_minid} \\
      --mincov ${params.blast_custom_mincov} \\
      --threads ${params.threads} \\
      --out ${prefix}_${customDB.baseName}_\${blast_cmd}.txt \\
      > ${prefix}_${customDB.baseName}_\${blast_cmd}.summary.txt ;

  # Step 4 - Get BED from blast
  awk \\
      '{print \$1 "\t" \$2 "\t" \$3}' \\
      ${prefix}_${customDB.baseName}_\${blast_cmd}.txt | \\
      tail -n +2 > \\
      ${prefix}_${customDB.baseName}_\${blast_cmd}.bed ;

  # Step 4 - Find intersection with annotation
  bedtools \\
      intersect \\
      -wa \\
      -a $gff \\
      -b ${prefix}_${customDB.baseName}_\${blast_cmd}.bed \\
      > ${prefix}_${customDB.baseName}_\${blast_cmd}.gff ;
  """
}
