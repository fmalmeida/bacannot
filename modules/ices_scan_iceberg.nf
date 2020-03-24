process iceberg {
  publishDir "${params.outdir}/${prefix}/ICEs", mode: 'copy'
  container = 'fmalmeida/bacannot:latest'
  tag "Scanning for ICEs with ICEberg database"

  input:
  tuple val(prefix), file(genes)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("blast_ices_iceberg.tsv")

  script:
  """
  # Prediction using the genes predicted from prokka

  diamond blastx --query-gencode 11 --db /work/iceberg/iceberg_prot -o blastprot2_result.tmp \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore stitle \
  --query $genes --query-cover ${params.diamond_MGEs_queryCoverage} ;

  ## Group it all
  awk -v id=${params.diamond_MGEs_identity} '{if (\$3 >= id) print }' \
  blastprot2_result.tmp > blast_ices_iceberg.tsv
  """
}
