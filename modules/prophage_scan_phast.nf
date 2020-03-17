process phast {
  publishDir "${params.outdir}/${prefix}/prophages", mode: 'copy'
  container = 'fmalmeida/bacannot:latest'
  tag "Scanning prophages with PHAST database"

  input:
  file genes
  val(prefix)

  output:
  file "blast_prophage_phast.tsv"

  script:
  """
  # Prediction using the genes predicted from prokka

  diamond blastx --query-gencode 11 --db /work/phast/phast_prot -o blast_result_genes.tmp \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore stitle \
  --query $genes --query-cover ${params.diamond_MGEs_queryCoverage} ;

  # Filter by identity

  awk -v id=${params.diamond_MGEs_identity} '{if (\$3 >= id) print }' \
  blast_result_genes.tmp > blast_prophage_phast.tsv
  """
}
