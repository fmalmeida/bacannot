process KEGG_DECODER {
  publishDir "${params.output}/${prefix}/KOfamscan", mode: 'copy'
  tag "${prefix}"
  label = [ 'misc', 'process_low' ]

  input:
  tuple val(prefix), path('input_mapper.txt')

  output:
  // Grab all outputs
  path("*") // Get all files to input directory
  tuple val(prefix), path("*.svg") // get svg

  script:
  """
  # Activate env
  export PATH=/opt/conda/envs/KEGGDecoder/bin:\$PATH

  # draw static heatmap
  KEGG-decoder \\
      --input input_mapper.txt \\
      --output ${prefix}_kegg-decoder_heatmap_static.tsv \\
      --vizoption static ;
  """
}
