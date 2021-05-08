process draw_GIs {
  publishDir "${params.outdir}/${prefix}/genomic_islands", mode: 'copy', saveAs: { filename ->
    if (filename == "plots") "$filename"
    else null
  }
  tag "Drawing predicted genomic islands"
  label 'main'

  input:
  tuple val(prefix), file(gff), file(gis_bed)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("plots") optional true
  tuple val(prefix), file("teste.png") optional true

  script:
  """
  # Create output directories
  mkdir -p plots plots/id_label plots/product_label ;

  # Get required files
  cp /work/bscripts/draw_gis.sh . ;
  cp /work/bscripts/input.fofn . ;

  # Draw genomic islands
  ./draw_gis.sh -i $gis_bed -g $gff -f input.fofn ;

  # Get one image
  name=\$(ls plots/product_label | head -n 1)
  [[ \$(ls plots/product_label/) ]] && cp "plots/product_label/\${name}" ./teste.png || echo "empty" ;
  """
}
