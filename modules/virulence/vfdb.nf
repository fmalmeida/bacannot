process VFDB {
  publishDir "${params.output}/${prefix}/virulence/vfdb", mode: 'copy'
  tag "${prefix}"
  label = [ 'misc', 'process_low' ]

  input:
  tuple val(prefix), file(genes)
  file(bacannot_db)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), path("${prefix}_vfdb_blastn_onGenes.summary.txt")
  tuple val(prefix), path("${prefix}_vfdb_blastn_onGenes.txt")
  path('*.txt')

  script:
  """
  # VFDB has nucleotide database
  run_blasts.py \\
      blastn \\
      --query $genes \\
      --db ${bacannot_db}/vfdb_db/sequences \\
      --minid ${params.blast_virulence_minid} \\
      --mincov ${params.blast_virulence_mincov} \\
      --threads ${params.threads} \\
      --out ${prefix}_vfdb_blastn_onGenes.txt \\
      --2way | \\
  sed -e 's/ACCESSION/VFDB_ID/g' > ${prefix}_vfdb_blastn_onGenes.summary.txt ;
  """
}
