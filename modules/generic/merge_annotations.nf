process merge_annotations {
  publishDir "${params.outdir}/${prefix}/gffs", mode: 'copy'
  label 'renv'
  tag "Merging all the different annotations"

  input:
  tuple val(prefix), file(draft), file("prokka_gff"), file(mlst), file(barrnap),
        file(gc_bedGraph), file(gc_chrSizes), file(kofamscan), file(vfdb),
        file(victors), file(amrfinder), file(rgi), file(iceberg), file(phast),
        file(phigaro), file(genomic_islands), file("digIS.gff")

  output:
  tuple val(prefix), file("${prefix}.gff")                  // Get main gff file
  file "virulence_vfdb.gff"           optional true         // Get VFDB virulence file
  file "ices_iceberg.gff"             optional true         // Get ICEberg ices file
  file "prophages_phast.gff"          optional true         // Get PHAST prophages file
  file "resistance_card.gff"          optional true         // Get CARD resistance file
  file "resistance_amrfinderplus.gff" optional true         // Get NDARO resistance file
  file "*.gff"                        optional true         // Get all subsets
  file "virulence_victors.gff"        optional true         // Get Victors virulence file

  script:
  """
  # Rename gff and remove sequence entries
  grep "ID=" prokka_gff > ${prefix}.gff ;

  ## Increment GFF with custom annotations
  ### VFDB
  [ ! -s ${prefix}_vfdb_blastn_onGenes.txt ] || addBlast2Gff.R -i $vfdb -g ${prefix}.gff -o ${prefix}.gff -d VFDB -t Virulence && \
  [ \$(grep "VFDB" ${prefix}.gff | wc -l) -eq 0 ] || grep "VFDB" ${prefix}.gff > virulence_vfdb.gff ;

  ### Victors
  [ ! -s ${prefix}_victors_blastp_onGenes.txt ] || addBlast2Gff.R -i $victors -g ${prefix}.gff -o ${prefix}.gff -d Victors -t Virulence && \
  [ \$(grep "Victors" ${prefix}.gff | wc -l) -eq 0 ] || grep "Victors" ${prefix}.gff > virulence_victors.gff ;

  ### KEGG Orthology
  ## Reformat KOfamscan Output
  [ ! -s ${prefix}_ko_forKEGGMapper.txt ] || awk -F'\\t' -v OFS='\\t' '{x=\$1;\$1="";a[x]=a[x]\$0}END{for(x in a)print x,a[x]}' $kofamscan  | \
                           sed -e 's/\\t/,/g' -e 's/,,/\\t/g' | awk  '\$2!=""' > formated.txt ;
  [ ! -s formated.txt ] || addKO2Gff.R -i formated.txt -g ${prefix}.gff -o ${prefix}.gff -d KEGG ;

  ### ICEs
  [ ! -s ${prefix}_iceberg_blastp_onGenes.txt ] || addBlast2Gff.R -i $iceberg -g ${prefix}.gff -o ${prefix}.gff -d ICEberg -t ICE && \
  [ \$(grep "ICEberg" ${prefix}.gff | wc -l) -eq 0 ] || grep "ICEberg" ${prefix}.gff > ices_iceberg.gff ;

  ### Prophages
  [ ! -s ${prefix}_phast_blastp_onGenes.txt ] || addBlast2Gff.R -i $phast -g ${prefix}.gff -o ${prefix}.gff -d PHAST -t Prophage && \
  [ \$(grep "PHAST" ${prefix}.gff | wc -l) -eq 0 ] || grep "PHAST" ${prefix}.gff > prophages_phast.gff ;

  ### Resistance
  #### RGI
  [ \$(cat RGI_${prefix}.txt | wc -l) -eq 1 ] || addRGI2gff.R -g ${prefix}.gff -i $rgi -o ${prefix}.gff ;
  [ \$(grep "CARD" ${prefix}.gff | wc -l) -eq 0 ] || grep "CARD" ${prefix}.gff > resistance_card.gff ;

  #### AMRFinderPlus
  [ \$(cat AMRFinder_resistance-only.tsv | wc -l) -eq 1 ] || addNCBIamr2Gff.R -g ${prefix}.gff -i $amrfinder -o ${prefix}.gff -t Resistance -d AMRFinderPlus ;
  [ \$(grep "AMRFinderPlus" ${prefix}.gff | wc -l) -eq 0 ] || grep "AMRFinderPlus" ${prefix}.gff > resistance_amrfinderplus.gff ;

  ### digIS transposable elements
  [ ! -s digIS.gff ] || cat ${prefix}.gff | sed 's/id=/ID=/g' | bedtools sort > tmp.gff && cat tmp.gff > ${prefix}.gff && rm tmp.gff ;
  """
}
