process refseq_masher {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "refseq_masher/$filename"
  }
  tag "${prefix}"
  label 'main'

  input:
  tuple val(prefix), file(genome)

  output:
  // Grab results
  tuple val(prefix), file("refseq_masher_results.txt")
  file("*_version.txt")

  script:
  """
  # activate env
  source activate PY36_env ;

  # Get tool version
  refseq_masher --version > refseq_masher_version.txt ;

  # Run tool
  refseq_masher -vv matches --top-n-results 10 --output-type tab $genome > refseq_masher_results.txt
  """
}
