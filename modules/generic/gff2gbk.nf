process gff2gbk {
  publishDir "${params.outdir}/${prefix}/gbk", mode: 'copy'
  tag "${prefix}"

  input:
  tuple val(prefix), file(gff), file(input)

  output:
  file "*.genbank"

  """
  seqret -sequence $input -feature -fformat gff -fopenfile $gff -osformat genbank \
  -osname_outseq ${prefix} -ofdirectory_outseq gbk_file -auto
  """
}
