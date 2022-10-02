process GFF2GBK {
  publishDir "${params.output}/${prefix}/gbk", mode: 'copy'
  tag "${prefix}"
  label = [ 'misc', 'process_ultralow' ]

  input:
  tuple val(prefix), file(gff), file(input)

  output:
  path "*.genbank", emit: results

  """
  # Activate env
  export PATH=/opt/conda/envs/antismash/bin:\$PATH

  # Run emboss seqret
  seqret \\
    -sequence $input \\
    -feature \\
    -fformat gff \\
    -fopenfile $gff \\
    -osformat genbank \\
    -osname_outseq ${prefix} \\
    -ofdirectory_outseq gbk_file \\
    -auto
  """
}
