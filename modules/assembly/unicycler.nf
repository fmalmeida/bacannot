process UNICYCLER {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "unicycler_${prefix}") "assembly"
    else null
  }
  label = [ 'process_high', 'error_retry' ]
  tag "${prefix}"

  conda "bioconda::unicycler=0.4.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/unicycler:0.4.8--py38h8162308_3' :
        'quay.io/biocontainers/unicycler:0.4.8--py38h8162308_3' }"

  input:
  tuple val(prefix), val(entrypoint), file(sread1), file(sread2), file(sreads), file(lreads), val(lr_type), file(fast5), val(assembly), val(resfinder_species)

  output:
  path "unicycler_${prefix}", emit: all // Save everything
  // Keep tuple structure to mixing channels
  tuple val("${prefix}"), val("${entrypoint}"), val("${sread1}"), val("${sread2}"), val("${sreads}"), path("${lreads}"), val("${lr_type}"), path("${fast5}"), path("unicycler_${prefix}.fasta"), val("${resfinder_species}"), emit: results
  path('unicycler_version.txt'), emit: version
  
  script:
  unpaired_param = (sreads.getName() != "input.3") ? "-s $sreads" : ""
  paired_param = (sread1.getName() != "input.1" && sread2.getName() != "input.2") ? "-1 $sread1 -2 $sread2" : ""
  lr_param = (lreads.getName() != "input.4") ? "-l $lreads" : ""

  """
  # Save unicycler version
  unicycler --version > unicycler_version.txt

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
