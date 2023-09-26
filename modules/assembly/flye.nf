process FLYE {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "flye_${prefix}") "assembly"
    else null
  }
  label 'process_high'
  tag "${prefix}"

  input:
  tuple val(prefix), val(entrypoint), file(sread1), file(sread2), file(sreads), file(lreads), val(lr_type), file(fast5), val(assembly), val(resfinder_species)

  output:
  path "flye_${prefix}", emit: all // Saves all files
  // Keep tuple structure to mixing channels
  tuple val("${prefix}"), val("${entrypoint}"), val("${sread1}"), val("${sread2}"), val("${sreads}"), path("${lreads}"), val("${lr_type}"), path("${fast5}"), path("flye_${prefix}.fasta"), val("${resfinder_species}"), emit: results
  path('flye_version.txt'), emit: version

  script:
  lr = (lr_type == 'nanopore') ? '--nano-raw' : '--pacbio-raw'
  dedup_lr = params.enable_deduplication ? 
    "gunzip -cf $lreads | awk '{if(NR%4==1) \$0=sprintf(\"@1_%d\",(1+i++)); print;}' | gzip -c > ${prefix}_deduplicated_reads.fastq.gz" :
    "ln -s $lreads ${prefix}_deduplicated_reads.fastq.gz"
  
  """
  # Save flye version
  flye -v > flye_version.txt ;

  # remove duplicate reads
  $dedup_lr

  # Run flye
  flye \\
    ${lr} \\
    ${prefix}_deduplicated_reads.fastq.gz \\
    --out-dir flye_${prefix} \\
    --threads $task.cpus &> flye.log ;

  # Save a copy for annotation
  cp flye_${prefix}/assembly.fasta flye_${prefix}.fasta
  """
}
