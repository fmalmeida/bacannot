process phigaro {
  publishDir "${params.outdir}/${prefix}/prophages/phigaro", mode: 'copy'
  tag "Scanning putative prophage sequences with phigaro"
  label 'main'

  input:
  tuple val(prefix), file("assembly.fasta")

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_phigaro.tsv")
  tuple val(prefix), file("${prefix}_phigaro.bed")
  tuple val(prefix), file("${prefix}_phigaro.html")

  script:
  """
  # Run phigaro
  phigaro -f assembly.fasta --config /work/phigaro_config.yml -t ${params.threads} \
  -e html tsv -o out.phg --delete-shorts -p --not-open ;

  # move
  mv out.phg/assembly.phigaro.tsv ${prefix}_phigaro.tsv ;
  mv out.phg/assembly.phigaro.html ${prefix}_phigaro.html;

  # Create BED
  grep -v "taxonomy" ${prefix}_phigaro.tsv | \
  awk 'BEGIN { FS = "\t"; OFS="\\t" } { print \$1,\$2,\$3 }' > ${prefix}_phigaro.bed
  """
}
