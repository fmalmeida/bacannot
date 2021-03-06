process find_GIs {
  publishDir "${params.outdir}/${prefix}/genomic_islands", mode: 'copy'
  errorStrategy 'retry'
  maxRetries 5
  tag "Predicting Genomic Islands with IslandPath-DIMOB"
  label 'main'

  input:
  tuple val(prefix), file("annotation.gbk")

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("${prefix}_predicted_GIs.bed")

  script:
  """
  # Activate environment
  source activate find_GIs ;

  # Split genbank files
  python /usr/local/bin/splitgenbank.py annotation.gbk && rm annotation.gbk ;

  # Run islandpath in each
  for file in \$(ls *.gbk); do grep -q "CDS" \$file && Dimob.pl \$file \${file%%.gbk}_GIs.txt 2> dimob.err ; done

  # Aggregate them
  for GI in \$(ls *.txt); do \
    name="\${GI%%_GIs.txt}" ;
    awk -v contig=\$name 'BEGIN { FS = "\t"; OFS="\\t" } { print contig,\$2,\$3 }' \$GI >> ${prefix}_predicted_GIs.bed ; \
  done
  """
}
