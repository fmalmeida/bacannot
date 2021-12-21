process KEGG_DECODER {
  publishDir "${params.output}/${prefix}/KOfamscan", mode: 'copy'
  tag "${prefix}"
  label 'misc'

  input:
  tuple val(prefix), file('input_mapper.txt')

  output:
  // Grab all outputs
  file("*") // Get all files to input directory
  tuple val(prefix), file("*.svg") // get svg

  script:
  """
  # draw static heatmap
  KEGG-decoder \\
      --input input_mapper.txt \\
      --output ${prefix}_kegg-decoder_heatmap_static.tsv \\
      --vizoption static ;
  """
}
