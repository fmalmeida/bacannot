process gff2gbk {
  publishDir "${params.outdir}/${prefix}/genbankFile", mode: 'copy'
  container = 'fmalmeida/bacannot:latest'

  input:
  tuple val(prefix), file(gff), file(input)

  output:
  file "*.genbank"

  """
  seqret -sequence $input -feature -fformat gff -fopenfile $gff -osformat genbank \
  -osname_outseq ${prefix} -ofdirectory_outseq gbk_file -auto
  """
}
