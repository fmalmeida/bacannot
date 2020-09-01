process find_GIs {
  publishDir "${params.outdir}/${prefix}/genomic_islands", mode: 'copy'
  tag "Predicting Genomic Islands with IslandPath-DIMOB"
  label 'main'

  input:
  tuple val(prefix), file("annotation.gbk")

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_predicted_GIs.bed")

  script:
  """
  source activate find_GIs ;
  python /work/pythonScripts/splitgenbank.py annotation.gbk && rm annotation.gbk ;
  for file in \$(ls *.gbk); do grep -q "CDS" \$file && Dimob.pl \$file \${file%%.gbk}_GIs.txt 2> dimob.err ; done
  for GI in \$(ls *.txt); do \
    name="\${GI%%_GIs.txt}" ;
    awk -v contig=\$name 'BEGIN { FS = "\t"; OFS="\\t" } { print contig,\$2,\$3 }' \$GI >> ${prefix}_predicted_GIs.bed ; \
  done
  """
}
