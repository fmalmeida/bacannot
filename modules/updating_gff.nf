process update_gff {
  publishDir "${params.outdir}/${prefix}/gffs", mode: 'copy'
  container = 'fmalmeida/bacannot:renv'

  input:
  tuple val(prefix), file(draft), file(gff), file(mlst), file(barrnap),
        file(gc_bedGraph), file(gc_chrSizes), file(kofamscan), file(vfdb),
        file(victors), file(amrfinder), file(rgi), file(iceberg), file(phast),
        file(phigaro), file(genomic_islands)

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
  mv $gff ${prefix}.gff ;


  ## Reformat KOfamscan Output
  while read line ; do id=\$(echo \$line | awk '{print \$1}') ; ko=\$(echo \$line | awk '{\$1=""; print \$0}' | \
  sed 's/\\s//' | sed 's/\\s/,/g') ; echo -e "\$id\\t\$ko" ; done < $kofamscan > formated.txt ;

  ## Increment GFF with custom annotations
  ### VFDB
  [ ! -s ${vfdb} ] || addBlast2Gff.R -i $vfdb -g ${prefix}.gff -o ${prefix}.gff -d VFDB -t Virulence && \
                            grep "VFDB" ${prefix}.gff > virulence_vfdb.gff ;

  ### Victors
  [ ! -s ${victors} ] || addBlast2Gff.R -i $victors -g ${prefix}.gff -o ${prefix}.gff -d Victors -t Virulence && \
                            grep "Victors" ${prefix}.gff > virulence_victors.gff ;

  ### KEGG Orthology
  [ ! -s ${kofamscan} ] || addKO2Gff.R -i formated.txt -g ${prefix}.gff -o ${prefix}.gff -d KEGG ;

  ### ICEs
  [ ! -s ${iceberg} ] || addBlast2Gff.R -i $iceberg -g ${prefix}.gff -o ${prefix}.gff -d ICEberg -t ICE && \
                               grep "ICEberg" ${prefix}.gff > ices_iceberg.gff ;

  ### Prophages
  [ ! -s ${phast} ] || addBlast2Gff.R -i $phast -g ${prefix}.gff -o ${prefix}.gff -d PHAST -t Prophage && \
                             grep "PHAST" ${prefix}.gff > prophages_phast.gff ;

  ### Resistance
  #### RGI
  [ ! -s ${rgi} ] || addRGI2gff.R -g ${prefix}.gff -i $rgi -o ${prefix}.gff && \
                                grep "CARD" ${prefix}.gff > resistance_card.gff ;

  #### AMRFinderPlus
  [ ! -s ${amrfinder} ] || addNCBIamr2Gff.R -g ${prefix}.gff -i $amrfinder -o ${prefix}.gff -t Resistance -d AMRFinderPlus && \
                                      grep "AMRFinderPlus" ${prefix}.gff > resistance_amrfinderplus.gff ;
  """
}
