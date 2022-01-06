process FIND_GIS {
  publishDir "${params.output}/${prefix}/genomic_islands", mode: 'copy'
  errorStrategy 'retry'
  maxRetries 5
  tag "${prefix}"
  label 'main'

  input:
  tuple val(prefix), file("annotation.gbk")

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_predicted_GIs.bed")

  script:
  """
  # activate env
  source activate PERL_env ;

  # Split genbank files
  splitgenbank.py annotation.gbk && rm annotation.gbk ;

  # Run islandpath in each
  touch ${prefix}_predicted_GIs.bed ;
  for file in \$(ls *.gbk); do \
    touch \${file%%.gbk}_GIs.txt ;
    grep -q "CDS" \$file && Dimob.pl \$file \${file%%.gbk}_GIs.txt 2> dimob.err ;
    name="\${file%%.gbk}" ;
    awk -v contig=\$name 'BEGIN { FS = "\\t"; OFS="\\t" } { print contig,\$2,\$3 }' \${file%%.gbk}_GIs.txt >> ${prefix}_predicted_GIs.bed ;
  done
  """
}
