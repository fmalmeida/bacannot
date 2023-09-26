process ISLANDPATH {
  publishDir "${params.output}/${prefix}/genomic_islands", mode: 'copy'
  tag "${prefix}"
  label = [ 'process_low' ]
  errorStrategy = 'retry'
  maxRetries    = 5

  input:
  tuple val(prefix), file("annotation.gbk")

  output:
  tuple val(prefix), path("${prefix}_predicted_GIs.bed"), emit: results

  script:
  """
  # Split genbank files
  splitgenbank.pl annotation.gbk && rm annotation.gbk ;

  # Run islandpath in each
  touch ${prefix}_predicted_GIs.bed ;
  for file in \$(ls *.gbk); do \
    touch \${file%%.gbk}_GIs.txt ;
    ( sed '/CDS.*::.*0/d' \$file | grep -q "CDS" ) && \\
        islandpath \\
        \$file \\
        \${file%%.gbk}_GIs.txt 2> dimob.err ;
    name="\${file%%.gbk}" ;
    awk -v contig=\$name 'BEGIN { FS = "\\t"; OFS="\\t" } { print contig,\$2,\$3 }' \${file%%.gbk}_GIs.txt >> ${prefix}_predicted_GIs.bed ;
  done
  """
}
