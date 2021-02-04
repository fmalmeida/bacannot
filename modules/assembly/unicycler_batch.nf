process unicycler_batch {
  publishDir "${params.outdir}/${id}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "unicycler_${id}") "assembly"
    else null
  }
  label 'assembly'
  tag { x }

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(sreads), file(lreads), val(lr_type), file(fast5), val(assembly), val(resfinder_species)

  output:
  file "unicycler_${id}" // Save everything
  // Keep tuple structure to mixing channels
  tuple val("${id}"), val("${entrypoint}"), val("${sread1}"), val("${sread2}"), val("${sreads}"), file("${lreads}"), val("${lr_type}"), file("${fast5}"), file("unicycler_${id}.fasta"), val("${resfinder_species}")
  file('unicycler_version.txt')

  when:
  if ((sread1.getName() != 'input.1' && sread2.getName() != 'input.2') || (sreads.getName() != 'input.3')) // Check if either paired or unpaired reads exists

  script:
  x = (lreads.getName() != "input.4") ? "Performing a hybrid assembly with Unicycler" : "Performing a illumina-only assembly with Unicycler"
  unpaired_param = (sreads.getName() != "input.3") ? "-s $sreads --no_correct" : ""
  paired_param = (sread1.getName() != "input.1" && sread2.getName() != "input.2") ? "-1 $sread1 -2 $sread2" : ""
  lr_param = (lreads.getName() != "input.4") ? "-l $lreads" : ""

  """
  # Save unicycler version
  unicycler --version > unicycler_version.txt

  # Run unicycler
  unicycler $paired_param $unpaired_param $lr_param -o unicycler_${id} -t ${params.threads} &> unicycler.log

  # Save copy for annotation
  cp unicycler_${id}/assembly.fasta unicycler_${id}.fasta
  """
}
