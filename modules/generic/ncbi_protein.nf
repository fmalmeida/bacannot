process NCBI_PROTEIN {
  publishDir "${params.output}/${prefix}/custom_annotations/NCBI_PROTEIN", mode: 'copy'
  tag "${prefix}"
  label 'main'

  input:
  tuple val(prefix), file(gff), file(proteins)
  file(ncbi_accs)

  output:
  // Grab results
  tuple val(prefix), val("NCBI_PROTEIN"), file("${prefix}_ncbi_protein_db_blastp.summary.txt"), file("${prefix}_ncbi_protein_db_blastp_blastp.gff")
  path("ncbi_protein_db.faa")
  path("${prefix}_ncbi_protein_db_blastp.txt")

  script:
  """
  # download and format ncbi protein entries for custom blastp
  for accession in \$(cat ${ncbi_accs}) ; do \\
    esearch -db protein -query "\${accession}" | \\
    efetch -format gp >> ncbi_protein_db.gbk; \\
  done

  # convert to formatted fasta
  gbk2faa.py ncbi_protein_db.gbk > ncbi_protein_db.faa

  # create blast db
  makeblastdb -in ncbi_protein_db.faa -dbtype prot -out customDB ;

  # execute blastp
  run_blasts.py \\
      blastp \\
      --query $proteins \\
      --db customDB \\
      --minid ${params.blast_custom_minid} \\
      --mincov ${params.blast_custom_mincov} \\
      --threads ${params.threads} \\
      --out ${prefix}_ncbi_protein_db_blastp.txt > \\
      ${prefix}_ncbi_protein_db_blastp.summary.txt ;

  # get BED from blastp
  awk \\
      '{print \$1 "\t" \$2 "\t" \$3}' \\
      ${prefix}_ncbi_protein_db_blastp_blastp.txt | \\
      tail -n +2 > ${prefix}_ncbi_protein_db_blastp_blastp.bed ;

  # find intersection with annotation
  bedtools \\
      intersect \\
      -wa \\
      -a $gff \\
      -b ${prefix}_ncbi_protein_db_blastp_blastp.bed > \\
      ${prefix}_ncbi_protein_db_blastp_blastp.gff ;
  """
}
