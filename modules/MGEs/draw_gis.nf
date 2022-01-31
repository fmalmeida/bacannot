process DRAW_GIS {
  publishDir "${params.output}/${prefix}/genomic_islands", mode: 'copy', saveAs: { filename ->
    if (filename == "plots") "$filename"
    else null
  }
  tag "${prefix}"
  

  input:
  tuple val(prefix), file(gff), file(gis_bed)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("plots")     optional true
  tuple val(prefix), file("teste.png") optional true

  script:
  """
  # create output directories
  mkdir -p plots plots/id_label plots/product_label ;

  # get required files
  cp /work/bscripts/draw_gis.sh . ;
  cp /work/bscripts/input.fofn . ;

  # draw genomic islands
  ./draw_gis.sh -i $gis_bed -g $gff -f input.fofn ;

  # get one image
  name=\$(ls plots/product_label | head -n 1)
  [[ \$(ls plots/product_label/) ]] && cp "plots/product_label/\${name}" ./teste.png || echo "empty" ;
  """
}
