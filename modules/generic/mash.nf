process REFSEQ_MASHER {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "refseq_masher/$filename"
  }
  tag "${prefix}"
  label = [ 'python', 'process_low' ]

  input:
  tuple val(prefix), path(genome)

  output:
  // Grab results
  tuple val(prefix), path("refseq_masher_results.txt")
  path("*_version.txt")

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
