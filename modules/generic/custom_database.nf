process CUSTOM_DATABASE {
  publishDir "${params.output}/${prefix}/custom_annotations/${customDB.baseName}", mode: 'copy'
  tag "${prefix}"
  label 'main'

  input:
  tuple val(prefix), file(gff), file(genome)
  each file(customDB)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), val("${customDB.baseName}"), path("${prefix}_${customDB.baseName}*.summary.txt"), path("${prefix}_${customDB.baseName}*.gff")
  path('*.txt') // Grab all
  path(customDB)

  script:
  """
  # Step 1 - Check if input is nucl or protein and prepare db
  if [ \$(grep -i "Protein" <(seqkit stats ${customDB}) | wc -l) -gt 0 ]
  then
        export blast_cmd="tblastn" ;
  else
        export blast_cmd="blastn" ;
  fi

  # Step 2 - Execute blast
  run_blasts.py \\
      \${blast_cmd} \\
      --query ${genome} \\
      --db ${customDB} \\
      --minid ${params.blast_custom_minid} \\
      --mincov ${params.blast_custom_mincov} \\
      --threads ${params.threads} \\
      --out ${prefix}_${customDB.baseName}_\${blast_cmd}.txt \\
      > ${prefix}_${customDB.baseName}_\${blast_cmd}.summary.txt ;

  # Step 3 - Get BED from blast
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
