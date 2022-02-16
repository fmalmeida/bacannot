process PLATON {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "platon") "plasmids/$filename"
    else null
  }
  tag "${prefix}"
  label 'python'

  input:
  tuple val(prefix), file(genome)
  file(bacannot_db)

  output:
  path("platon")
  tuple val(prefix), path("platon/${prefix}.tsv")
  path("platon_version.txt")

  script:
  """
  # Get version
  platon --version > platon_version.txt ;

  # Run platon
  platon \\
      --db ${bacannot_db}/platon_db/ \\
      --output platon \\
      --threads ${params.threads} \\
      $genome > tmp.txt || true ;
  [ -s platon/${prefix}.tsv ] || cat tmp.txt > platon/${prefix}.tsv ;
  """
}
