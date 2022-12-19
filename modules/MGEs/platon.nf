process PLATON {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "platon") "plasmids/$filename"
    else null
  }
  tag "${prefix}"
  label = [ 'python', 'process_medium' ]

  input:
  tuple val(prefix), file(genome)
  file(bacannot_db)

  output:
  tuple val(prefix), path("platon")              , emit: all
  tuple val(prefix), path("platon/${prefix}.tsv"), emit: results
  path("platon_version.txt")                     , emit: version

  script:
  """
  # Get version
  platon --version > platon_version.txt ;

  # Run platon
  platon \\
      --db ${bacannot_db}/platon_db/ \\
      --output platon \\
      --threads $task.cpus \\
      $genome > tmp.txt || true ;
  [ -s platon/${prefix}.tsv ] || cat tmp.txt > platon/${prefix}.tsv ;
  """
}
