process MERGE_ANNOTATIONS {
  publishDir "${params.output}/${prefix}/gffs", mode: 'copy'
  label = [ 'renv', 'process_medium', 'error_retry' ]
  tag "${prefix}"

  input:
  tuple val(prefix), file('prokka_gff'), file(kofamscan), file(vfdb), file(victors), file(amrfinder), file(resfinder), file(rgi), file(iceberg), file(phast), file('digis_gff'), file(custom_databases)

  output:
  tuple val(prefix), path("${prefix}.gff")                      , emit: gff
  tuple val(prefix), path("transposable_elements_digis.gff")    , emit: digis_gff
  tuple val(prefix), path("custom_database_*.gff") optional true, emit: customdb_gff
  path("*.gff")                                                 , emit: all

  script:
  """
  # Rename gff and remove sequence entries
  # bakta has region entries
  awk '\$3 == "CDS"' prokka_gff | grep "ID=" > ${prefix}.gff ;

  ## Increment GFF with custom annotations
  ### VFDB
  if [ ! \$(cat $vfdb | wc -l) -le 1 ]
  then
    addBlast2Gff.R -i $vfdb -g ${prefix}.gff -o ${prefix}.gff -d VFDB -t Virulence ;
    grep "VFDB" ${prefix}.gff > virulence_vfdb.gff ;
  fi

  ### Victors
  if [ ! \$(cat $victors | wc -l) -le 1 ]
  then 
    addBlast2Gff.R -i $victors -g ${prefix}.gff -o ${prefix}.gff -d Victors -t Virulence ;
    grep "Victors" ${prefix}.gff > virulence_victors.gff ;
  fi

  ### KEGG Orthology
  ## Reformat KOfamscan Output
  if [ ! \$(cat $kofamscan | wc -l) -eq 0 ]
  then
    awk \\
      -F'\\t' \\
      -v OFS='\\t' \\
      '{x=\$1;\$1="";a[x]=a[x]\$0}END{for(x in a)print x,a[x]}' \\
      $kofamscan  | \\
    sed \\
      -e 's/\\t/,/g' \\
      -e 's/,,/\\t/g' | \\
    awk  '\$2!=""' > formated.txt ;
    addKO2Gff.R -i formated.txt -g ${prefix}.gff -o ${prefix}.gff -d KEGG ;
  fi

  ### ICEs
  if [ ! \$(cat $iceberg | wc -l) -le 1 ]
  then
    addBlast2Gff.R -i $iceberg -g ${prefix}.gff -o ${prefix}.gff -d ICEberg -t ICE ;
    grep "ICEberg" ${prefix}.gff > ices_iceberg.gff ;
  fi

  ### Prophages
  if [ ! \$(cat $phast | wc -l) -le 1 ]
  then
    addBlast2Gff.R -i $phast -g ${prefix}.gff -o ${prefix}.gff -d PHAST -t Prophage ;
    grep "PHAST" ${prefix}.gff > prophages_phast.gff ;
  fi

  ### Resistance
  #### RGI
  if [ ! \$(cat $rgi | wc -l) -le 1 ]
  then
    addRGI2gff.R -g ${prefix}.gff -i $rgi -o ${prefix}.gff ;
    grep "CARD" ${prefix}.gff > resistance_card.gff ;
  fi

  #### AMRFinderPlus
  if [ ! \$(cat $amrfinder | wc -l) -le 1 ]
  then 
    addNCBIamr2Gff.R -g ${prefix}.gff -i $amrfinder -o ${prefix}.gff -t Resistance -d AMRFinderPlus ;
    grep "AMRFinderPlus" ${prefix}.gff > resistance_amrfinderplus.gff ;
  fi

  #### Resfinder
  if [ ! \$(cat $resfinder | wc -l) -eq 0 ]
  then
    bedtools intersect -a $resfinder -b ${prefix}.gff -wo | sort -k19,19 -r | awk -F '\\t' '!seen[\$9]++ {print \$18}' > resfinder_intersected.txt ;
    addBedtoolsIntersect.R -g ${prefix}.gff -t resfinder_intersected.txt --type Resistance --source Resfinder -o ${prefix}.gff ;
    grep "Resfinder" ${prefix}.gff > resistance_resfinder.gff ;
    rm -f resfinder_intersected.txt ;
  fi

  #### Custom Blast databases
  for file in ${custom_databases.join(" ")} ;
  do
    if [ ! \$(cat \$file | wc -l) -eq 0 ]
    then
      db=\${file%%_custom_db.gff} ;
      bedtools intersect -a \${file} -b ${prefix}.gff -wo | sort -k19,19 -r | awk -F '\\t' '!seen[\$9]++ {print \$18}' > bedtools_intersected.txt ;
      addBedtoolsIntersect.R -g ${prefix}.gff -t bedtools_intersected.txt --type "CDS" --source "\${db}" -o ${prefix}.gff ;
      grep "\${db}" ${prefix}.gff > custom_database_\${db}.gff ;
      rm -f bedtools_intersected.txt ;
    fi
  done

  ### digIS transposable elements
  touch transposable_elements_digis.gff
  if [ -s digis_gff ]
  then
    ( cat digis_gff | sed 's/id=/ID=/g' > transposable_elements_digis.gff && rm digis_gff ) ;
    cat ${prefix}.gff transposable_elements_digis.gff | bedtools sort > tmp.out.gff ;
    ( cat tmp.out.gff > ${prefix}.gff && rm tmp.out.gff );
  fi
  """
}
