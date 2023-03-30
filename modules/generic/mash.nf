process REFSEQ_MASHER {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "refseq_masher/$filename"
  }
  tag "${prefix}"
  label = [ 'process_low' ]

  input:
  tuple val(prefix), path(genome)

  output:
  tuple val(prefix), path("refseq_masher_results.txt"), emit: results
  path("*_version.txt")                               , emit: version

  script:
  """
  # Get tool version
  refseq_masher --version > refseq_masher_version.txt ;

  # Run tool
  refseq_masher \\
    -vv matches \\
    --top-n-results 10 \\
    --output-type tab \\
    $genome > refseq_masher_results.txt
  """
}
