process unicycler {
  publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "unicycler_${prefix}") "assembly"
    else null
  }
  label 'assembly'
  tag "${prefix}"

  input:
  tuple val(id), file(sread1), file(sread2)
  file(sreads)
  file(lreads)

  output:
  file "unicycler_${prefix}/*"
  file("unicycler_${prefix}.fasta")
  file('unicycler_version.txt')

  script:
  unpaired_param = (sreads.getName() != "input.3") ? "-s $sreads" : ""
  paired_param = (sread1.getName() != "input.1" && sread2.getName() != "input.2") ? "-1 $sread1 -2 $sread2" : ""
  lr_param = (lreads.getName() != "input.4") ? "-l $lreads" : ""
  prefix = params.prefix
  """
  # Save unicycler version
  unicycler --version > unicycler_version.txt

  # Run unicycler
  unicycler $paired_param $unpaired_param $lr_param -o unicycler_${prefix} -t ${params.threads} &> unicycler.log

  # Save copy for annotation
  cp unicycler_${prefix}/assembly.fasta unicycler_${prefix}.fasta
  """
}
