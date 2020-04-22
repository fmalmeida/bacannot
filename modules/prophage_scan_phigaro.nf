process phigaro {
  publishDir "${params.outdir}/${prefix}/prophages/phigaro", mode: 'copy'
  container = 'fmalmeida/bacannot:latest'
  tag "Scanning putative prophages with phigaro"

  input:
  tuple val(prefix), file("assembly.fasta")

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_phigaro.phg")
  tuple val(prefix), file("${prefix}_phigaro.bed")
  tuple val(prefix), file("${prefix}_phigaro.phg.html")

  script:
  """
  touch ${prefix}_phigaro.phg ${prefix}_phigaro.phg.html ;

  # Run phigaro
  phigaro -f assembly.fasta --config /work/phigaro_config.yml \
  -e html txt -o ${prefix}_phigaro.phg --delete-shorts -p --not-open ;

  # Create BED
  grep -v "taxonomy" ${prefix}_phigaro.phg | \
  awk 'BEGIN { FS = "\t"; OFS="\\t" } { print \$1,\$2,\$3 }' > ${prefix}_phigaro.bed
  """
}
