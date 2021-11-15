process platon {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "platon") "plasmids/$filename"
    else null
  }
  tag "${prefix}"
  label 'main'

  input:
  tuple val(prefix), file(genome)

  output:
  file("platon")
  tuple val(prefix), file("platon/${prefix}.tsv")
  file("platon_version.txt")

  script:
  """
  # Get version
  platon --version > platon_version.txt ;

  # Unpack database
  tar zxvf /work/platon/db.tar.gz ;

  # Run platon
  platon --db db/ --output platon --threads ${params.threads} $genome > tmp.txt || true ;
  [ -s platon/${prefix}.tsv ] || cat tmp.txt > platon/${prefix}.tsv ;
  """
}
