process flye {
  publishDir "${params.outdir}/${id}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "flye_${id}") "assembly"
    else null
  }
  label 'assembly'
  tag "${id}"

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(sreads), file(lreads), val(lr_type), file(fast5), val(assembly), val(resfinder_species)

  output:
  file "flye_${id}" // Saves all files
  // Keep tuple structure to mixing channels
  tuple val("${id}"), val("${entrypoint}"), val("${sread1}"), val("${sread2}"), val("${sreads}"), file("${lreads}"), val("${lr_type}"), file("${fast5}"), file("flye_${id}.fasta"), val("${resfinder_species}")
  file('flye_version.txt')

  script:
  lr = (lr_type == 'nanopore') ? '--nano-raw' : '--pacbio-raw'
  """
  source activate flye ;

  # Save flye version
  flye -v > flye_version.txt ;

  # Run flye
  flye ${lr} $lreads --plasmids --out-dir flye_${id} --threads ${params.threads} &> flye.log ;

  # Save a copy for annotation
  cp flye_${id}/assembly.fasta flye_${id}.fasta
  """
}
