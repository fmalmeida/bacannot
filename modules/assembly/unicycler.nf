process unicycler {
  publishDir "${params.outdir}/assembly", mode: 'copy'
  label 'assembly'
  tag { x }

  input:
  tuple val(id), file(sread1), file(sread2)
  file(sreads)
  file(lreads)

  output:
  file "unicycler_${id}" // Save everything
  file("unicycler_${id}.fasta")

  script:
  x = "Performing a illumina-only assembly with Unicycler"
  if ((params.sreads_single) && (params.sreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sreads --no_correct"
  } else if ((params.sreads_single) && (!params.sreads_paired)) {
    parameter = "-s $sreads --no_correct"
    id = sreads.getSimpleName()
  } else if ((params.sreads_paired) && (!params.sreads_single)) {
    parameter = "-1 $sread1 -2 $sread2"
  }

  if (params.lreads) {
    lr_param = "-l $lreads"
    x = "Performing a hybrid assembly with Unicycler"
  } else {
    lr_param = ""
  }

  """
  unicycler $parameter $lr_param -o unicycler_${id} -t ${params.threads} &> unicycler.log

  # get copy
  cp unicycler_${id}/assembly.fasta unicycler_${id}.fasta
  """
}
