process kegg_decoder {
  publishDir "${params.output}/${prefix}/KOfamscan", mode: 'copy'
  tag "${prefix}"
  label 'kofam'

  input:
  tuple val(prefix), file('input_mapper.txt')

  output:
  // Grab all outputs
  file("*") // Get all files to input directory
  tuple val(prefix), file("*.svg") // get svg

  script:
  """
  # KEGG-DECODER
  source activate kegg-decoder-env ;

  # Draw static heatmap
  KEGG-decoder --input input_mapper.txt\
  --output ${prefix}_kegg-decoder_heatmap-static.tsv --vizoption static ;
  """
}
