process find_GIs {
  publishDir "${params.outdir}/${prefix}/predicted_GIs", mode: 'copy'
  container = 'fmalmeida/bacannot:latest'
  tag "Predicting Genomic Islands with IslandPath-DIMOB"

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
    awk -v contig="\$( echo \"gnl|${params.prokka_center}|\${GI%%_GIs.txt}\" )" \
    'BEGIN { FS = "\t"; OFS="\\t" } { print contig,\$2,\$3 }' \$GI >> ${prefix}_predicted_GIs.bed ; \
  done
  """
}
