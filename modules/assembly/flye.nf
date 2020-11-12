process flye {
  publishDir "${params.outdir}/assembly", mode: 'copy', overwrite: true
  label 'assembly'
  tag "Performing a longreads only assembly with Flye"

  input:
  file lreads

  output:
  file "flye_${lrID}" // Saves all files
  file("flye_${lrID}.fasta") // Gets contigs file

  script:
  lr = (params.lreads_type == 'nanopore') ? '--nano-raw' : '--pacbio-raw'
  lrID = lreads.getSimpleName()
  """
  source activate flye ;
  flye ${lr} $lreads --plasmids --out-dir flye_${lrID} --threads ${params.threads} &> flye.log ;
  cp flye_${lrID}/assembly.fasta flye_${lrID}.fasta
  """
}
