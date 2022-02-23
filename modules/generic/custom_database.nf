process CUSTOM_DATABASE {
  publishDir "${params.output}/${prefix}/custom_annotations/${customDB.baseName}", mode: 'copy'
  tag "${prefix}"
  label = [ 'misc', 'process_low' ]

  input:
  tuple val(prefix), file(gff), file(genome)
  each file(customDB)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), val("${customDB.baseName}"), path("${prefix}_${customDB.baseName}*.summary.txt")
  tuple val(prefix), path("${customDB.baseName}_custom_db.gff")
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

  # Step 3 - Produce custom gff
  tail -n+2 ${prefix}_${customDB.baseName}_\${blast_cmd}.summary.txt | \\
  awk \\
    -v source="${customDB.baseName}" \\
    -F'\\t' \\
    'BEGIN{ OFS="\\t"; }
    { 
        atts="Additional_database="\$10";"\$10"_Acc="\$11";"\$10"_Target="\$5";"\$10"_Product="\$12";"\$10"_Description="\$13;
        if (\$4 == "-") {
            print \$1,source,"CDS",\$3,\$2,".",\$4,"0",atts
        } else {
            print \$1,source,"CDS",\$2,\$3,".",\$4,"0",atts
        }
    }' | \\
  bedtools sort > ${customDB.baseName}_custom_db.gff
  """

}