process phigaro {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename == "out.phg") null
    else if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "prophages/phigaro/$filename"
  }
  tag "Scanning putative prophage sequences with phigaro"
  label 'main'

  input:
  tuple val(prefix), file("assembly.fasta")

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_phigaro.tsv")
  tuple val(prefix), file("${prefix}_phigaro.bed")
  tuple val(prefix), file("${prefix}_phigaro.html") optional true
  file('phigaro_version.txt')

  script:
  """
  # Get tool version
  phigaro -V > phigaro_version.txt ;

  # Run phigaro
  phigaro -f assembly.fasta --config /work/phigaro_config.yml -t ${params.threads} \
  -e html tsv -o out.phg --delete-shorts -p --not-open ;

  # Change names
  [ ! -s out.phg/assembly.phigaro.tsv  ] || mv out.phg/assembly.phigaro.tsv ${prefix}_phigaro.tsv ;
  [ ! -s out.phg/assembly.phigaro.html ] || mv out.phg/assembly.phigaro.html ${prefix}_phigaro.html ;

  # Create BED
  grep -v "taxonomy" ${prefix}_phigaro.tsv | \
  awk 'BEGIN { FS = "\t"; OFS="\\t" } { print \$1,\$2,\$3 }' > ${prefix}_phigaro.bed
  """
}
