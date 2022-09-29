process CUSTOM_DATABASE {
  publishDir "${params.output}/${prefix}/custom_annotations/${customDB.baseName}", mode: 'copy'
  tag "${prefix}"
  label = [ 'misc', 'process_low' ]

  input:
  tuple val(prefix), file(gff), file(genome)
  each file(customDB)

  output:
  tuple val(prefix), val("${customDB.baseName}"), path("${prefix}_${customDB.baseName}*.summary.txt"), emit: summary
  tuple val(prefix), path("${customDB.baseName}_custom_db.gff")                                      , emit: gff
  path('*.txt')                                                                                      , emit: all
  path(customDB)                                                                                     , emit: db

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
  if [[ \${blast_cmd} == "blastn" ]]
  then
      ## In case the blast dabase throw error below:
      ## BLAST Database error: No alias or index file found for nucleotide database
      makeblastdb \\
          -dbtype nucl \\
          -in ${customDB} \\
          -out ${customDB} ;
  fi

  if [[ \${blast_cmd} == "tblastn" ]]
  then
      ## In case the blast dabase throw error below:
      ## BLAST Database error: No alias or index file found for nucleotide database
      makeblastdb \\
          -dbtype prot \\
          -in ${customDB} \\
          -out ${customDB} ;
  fi

  run_blasts.py \\
      \${blast_cmd} \\
      --query ${genome} \\
      --db ${customDB} \\
      --minid ${params.blast_custom_minid} \\
      --mincov ${params.blast_custom_mincov} \\
      --threads $task.cpus \\
      --out ${prefix}_${customDB.baseName}_\${blast_cmd}_onGenome.txt \\
      > ${prefix}_${customDB.baseName}_\${blast_cmd}_onGenome.summary.txt ;

  # Step 3 - Produce custom gff
  tail -n+2 ${prefix}_${customDB.baseName}_\${blast_cmd}_onGenome.summary.txt | \\
  awk \\
    -v source="${customDB.baseName}" \\
    -F'\\t' \\
    'BEGIN{ OFS="\\t"; }
    { 
        atts="Additional_database="\$10";"\$10":Acc="\$11";"\$10":Target="\$5";"\$10":Product="\$12";"\$10":Description="\$13;
        if (\$4 == "-") {
            print \$1,source,"CDS",\$3,\$2,".",\$4,"0",atts
        } else {
            print \$1,source,"CDS",\$2,\$3,".",\$4,"0",atts
        }
    }' | \\
  bedtools sort > ${customDB.baseName}_custom_db.gff
  """

}
