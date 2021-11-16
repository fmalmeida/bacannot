process flye {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "flye_${prefix}") "assembly"
    else null
  }
  label 'assembly'
  tag "${prefix}"

  input:
  tuple val(prefix), val(entrypoint), file(sread1), file(sread2), file(sreads), file(lreads), val(lr_type), file(fast5), val(assembly), val(resfinder_species)

  output:
  file "flye_${prefix}" // Saves all files
  // Keep tuple structure to mixing channels
  tuple val("${prefix}"), val("${entrypoint}"), val("${sread1}"), val("${sread2}"), val("${sreads}"), file("${lreads}"), val("${lr_type}"), file("${fast5}"), file("flye_${prefix}.fasta"), val("${resfinder_species}")
  file('flye_version.txt')

  script:
  lr = (lr_type == 'nanopore') ? '--nano-raw' : '--pacbio-raw'
  """
  source activate flye ;

  # Save flye version
  flye -v > flye_version.txt ;

  # Run flye
  flye ${lr} $lreads --plasmids --out-dir flye_${prefix} --threads ${params.threads} &> flye.log ;

  # Save a copy for annotation
  cp flye_${prefix}/assembly.fasta flye_${prefix}.fasta
  """
}
