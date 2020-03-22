process update_gff {
  publishDir "${params.outdir}/${prefix}/gffs", mode: 'copy'
  container = 'fmalmeida/bacannot:renv'

  input:
  val(prefix)
  file 'gff'
  file 'kofamscan.txt'
  file 'vfdb_blast'
  file "AMRFinder_output.tsv"
  file 'RGI_output.txt'
  file 'iceberg_blast'
  file 'phast_blast'

  output:
  file "${prefix}.gff"                // Get main gff file
  file "virulence_vfdb.gff"           optional true         // Get VFDB virulence file
  file "ices_iceberg.gff"             optional true         // Get ICEberg ices file
  file "prophages_phast.gff"          optional true         // Get PHAST prophages file
  file "resistance_card.gff"          optional true         // Get CARD resistance file
  file "resistance_amrfinderplus.gff" optional true         // Get NDARO resistance file
  file "*.gff"                        // Get all subsets

  script:
  """
  # Rename gff file
  mv gff ${prefix}.gff ;


  ## Reformat KOfamscan Output
  while read line ; do id=\$(echo \$line | awk '{print \$1}') ; ko=\$(echo \$line | awk '{\$1=""; print \$0}' | \
  sed 's/\\s//' | sed 's/\\s/,/g') ; echo -e "\$id\\t\$ko" ; done < kofamscan.txt > formated.txt ;

  ## Increment GFF with custom annotations
  ### VFDB
  [ ! -s vfdb_blast ] || addBlast2Gff.R -i vfdb_blast -g ${prefix}.gff -o ${prefix}.gff -d VFDB -t Virulence && \
                            grep "VFDB" ${prefix}.gff > virulence_vfdb.gff ;

  ### KEGG Orthology
  [ ! -s kofamscan.txt ] || addKO2Gff.R -i formated.txt -g ${prefix}.gff -o ${prefix}.gff -d KEGG ;

  ### ICEs
  [ ! -s iceberg_blast ] || addBlast2Gff.R -i iceberg_blast -g ${prefix}.gff -o ${prefix}.gff -d ICEberg -t ICE && \
                               grep "ICEberg" ${prefix}.gff > ices_iceberg.gff ;

  ### Prophages
  [ ! -s phast_blast ] || addBlast2Gff.R -i phast_blast -g ${prefix}.gff -o ${prefix}.gff -d PHAST -t Prophage && \
                             grep "PHAST" ${prefix}.gff > prophages_phast.gff ;

  ### Resistance
  #### RGI
  [ ! -s RGI_output.txt ] || addRGI2gff.R -g ${prefix}.gff -i RGI_output.txt -o ${prefix}.gff && \
                                grep "CARD" ${prefix}.gff > resistance_card.gff ;

  #### AMRFinderPlus
  [ ! -s AMRFinder_output.tsv ] || addNCBIamr2Gff.R -g ${prefix}.gff -i AMRFinder_output.tsv -o ${prefix}.gff -t Resistance -d AMRFinderPlus && \
                                      grep "AMRFinderPlus" ${prefix}.gff > resistance_amrfinderplus.gff ;
  """
}
