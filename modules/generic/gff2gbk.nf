process GFF2GBK {
  publishDir "${params.output}/${prefix}/gbk", mode: 'copy'
  label 'main'
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
