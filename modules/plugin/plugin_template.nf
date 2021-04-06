process plugin_name {
  plugin_outdir = 'plugin_out'
  publishDir "${params.outdir}/${prefix}/${plugin_outdir}", mode: 'copy'
  tag "Drawing predicted genomic islands"
  //container ''

  input:
  tuple val(prefix), file(gff), file(gis_bed)

  output:

  script:
  """
  ls 
  """
}
