process flye {
  publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "flye_${prefix}") "assembly"
    else null
  }
  label 'assembly'
  tag "Performing a longreads only assembly with Flye"

  input:
  file(lreads)

  output:
  file "flye_${prefix}"
  file("flye_${prefix}.fasta")
  file('flye_version.txt')

  script:
  lr = (params.lreads_type == 'nanopore') ? '--nano-raw' : '--pacbio-raw'
  prefix = params.lreads_type
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
