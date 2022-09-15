process DRAW_GIS {
  publishDir "${params.output}/${prefix}/genomic_islands", mode: 'copy', saveAs: { filename ->
    if (filename == "plots") "$filename"
    else null
  }
  tag "${prefix}"
  label = [ 'misc', 'process_ultralow' ]
  

  input:
  tuple val(prefix), file(gff), file(gis_bed)

  output:
  tuple val(prefix), file("plots")     optional true, emit: all
  tuple val(prefix), file("teste.png") optional true, emit: example

  script:
  """
  # create output directories
  mkdir \\
    -p plots \\
    plots/id_label \\
    plots/product_label ;

  # draw genomic islands
  draw_gis.sh \\
    -i $gis_bed \\
    -g $gff \\
    -f \$(which input.fofn) ;

  # get one image
  name=\$(ls plots/product_label | head -n 1)
  [[ \$(ls plots/product_label/) ]] && \\
    cp "plots/product_label/\${name}" ./teste.png || \\
    echo "empty" ;
  """
}
