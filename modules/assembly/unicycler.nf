process UNICYCLER {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "unicycler_${prefix}") "assembly"
    else null
  }
  label 'process_high'
  tag "${prefix}"

  input:
  tuple val(prefix), val(entrypoint), file(sread1), file(sread2), file(sreads), file(lreads), val(lr_type), file(fast5), val(assembly), val(resfinder_species)

  output:
  path "unicycler_${prefix}", emit: all // Save everything
  // Keep tuple structure to mixing channels
  tuple val("${prefix}"), val("${entrypoint}"), val("${sread1}"), val("${sread2}"), val("${sreads}"), path("${lreads}"), val("${lr_type}"), path("${fast5}"), path("unicycler_${prefix}.fasta"), val("${resfinder_species}"), emit: results
  path('unicycler_version.txt'), emit: version
  
  script:
  unpaired_param = ""
  dedup_sreads   = ""
  paired_param   = ""
  dedup_paired   = ""
  lr_param       = ""
  dedup_lr       = ""

  // sreads
  if (sreads.getName() != "input.3") {

    dedup_sreads = params.enable_deduplication ? 
      "gunzip -cf $sreads | awk '{if(NR%4==1) \$0=sprintf(\"@1_%d\",(1+i++)); print;}' | gzip -c > ${prefix}_deduplicated_sreads.fastq.gz" :
      "ln -s $sreads ${prefix}_deduplicated_sreads.fastq.gz"
    
    unpaired_param = "-s ${prefix}_deduplicated_sreads.fastq.gz"

  }

  // paired
  if (sread1.getName() != "input.1" && sread2.getName() != "input.2") {

    dedup_paired = params.enable_deduplication ?  
      "gunzip -cf $sread1 | awk '{if(NR%4==1) \$0=sprintf(\"@1_%d\",(1+i++)); print;}' | gzip -c > ${prefix}_deduplicated_sread_R1.fastq.gz && gunzip -cf $sread2 | awk '{if(NR%4==1) \$0=sprintf(\"@1_%d\",(1+i++)); print;}' | gzip -c > ${prefix}_deduplicated_sread_R2.fastq.gz" :
      "ln -s $sread1 ${prefix}_deduplicated_sread_R1.fastq.gz && ln -s $sread2 ${prefix}_deduplicated_sread_R2.fastq.gz"
    
    paired_param   = "-1 ${prefix}_deduplicated_sread_R1.fastq.gz -2 ${prefix}_deduplicated_sread_R2.fastq.gz"

  }

  // lreads
  if (lreads.getName() != "input.4") {

    dedup_lr = params.enable_deduplication ? 
      "gunzip -cf $lreads | awk '{if(NR%4==1) \$0=sprintf(\"@1_%d\",(1+i++)); print;}' | gzip -c > ${prefix}_deduplicated_lreads.fastq.gz" :
      "ln -s $lreads ${prefix}_deduplicated_lreads.fastq.gz"
    
    lr_param   = "-l $lreads"

  }

  """
  # Save unicycler version
  unicycler --version > unicycler_version.txt

  # remove duplicate reads
  $dedup_sreads
  $dedup_paired
  $dedup_lr

  # Run unicycler
  unicycler \\
    $paired_param \\
    $unpaired_param \\
    $lr_param \\
    -o unicycler_${prefix} \\
    -t $task.cpus &> unicycler.log

  # Save copy for annotation
  cp unicycler_${prefix}/assembly.fasta unicycler_${prefix}.fasta
  """
}
