process vfdb {
  publishDir "${params.outdir}/${prefix}/virulence", mode: 'copy'
  container = 'fmalmeida/bacannot:latest'
  tag "Scanning virulence genes with VFDB"

  input:
  file genes
  val(prefix)

  output:
  file "blast_virulence_vfdb.tsv" // Save the results in tabular blast format

  script:
  """
  # Prediction using the genes predicted from prokka

  diamond blastx --query-gencode 11 --db /work/vfdb/vfdb_prot -o blast_result_genes.tmp \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore stitle \
  --query $genes --query-cover ${params.diamond_virulence_queryCoverage} ;

  # Filter by identity

  awk -v id=${params.diamond_virulence_identity} '{if (\$3 >= id) print }' \
  blast_result_genes.tmp > blast_virulence_vfdb.tsv
  """

}
