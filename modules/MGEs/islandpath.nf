process ISLANDPATH {
  publishDir "${params.output}/${prefix}/genomic_islands", mode: 'copy'
  tag "${prefix}"
  label 'perl'

  input:
  tuple val(prefix), file("annotation.gbk")

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), path("${prefix}_predicted_GIs.bed")

  script:
  """
  # Split genbank files
  splitgenbank.pl annotation.gbk && rm annotation.gbk ;

  # Run islandpath in each
  touch ${prefix}_predicted_GIs.bed ;
  for file in \$(ls *.gbk); do \
    touch \${file%%.gbk}_GIs.txt ;
    grep -q "CDS" \$file && \\
        islandpath \\
        \$file \\
        \${file%%.gbk}_GIs.txt 2> dimob.err ;
    name="\${file%%.gbk}" ;
    awk -v contig=\$name 'BEGIN { FS = "\\t"; OFS="\\t" } { print contig,\$2,\$3 }' \${file%%.gbk}_GIs.txt >> ${prefix}_predicted_GIs.bed ;
  done
  """
}
