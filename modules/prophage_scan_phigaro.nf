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
  touch assembly.phg assembly.phg.html ;

  # Filter fasta with minimum size of 20 kb (required by phigaro)
  seqtk seq -L 20000 assembly.fasta > assembly-L20000.fasta ;

  # Run phigaro
  phigaro -f assembly-L20000.fasta -c /work/phigaro/config.yml -e html txt -o ${prefix}_phigaro.phg -p --not-open ;

  # Create BED
  grep -v "taxonomy" ${prefix}_phigaro.phg | \
  awk 'BEGIN { FS = "\t"; OFS="\\t" } { print \$1,\$2,\$3 }' > ${prefix}_phigaro.bed
  """
}
