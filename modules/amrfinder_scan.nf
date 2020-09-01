process amrfinder {
  publishDir "${params.outdir}/${prefix}/resistance/AMRFinderPlus", mode: 'copy'
  tag "Scanning AMR genes with AMRFinderPlus"
  label 'main'

  input:
  tuple val(prefix), file(proteins)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("AMRFinder_resistance-only.tsv")
  tuple val(prefix), file("AMRFinder_complete.tsv")

  script:
  """
  source activate AMRFINDERPLUS ;
  amrfinder -p $proteins --plus -o AMRFinder_complete.tsv --threads ${params.threads} \
  --ident_min \$(echo "scale=2; ${params.diamond_resistance_minid}/100" | bc -l ) \
  --coverage_min \$(echo "scale=2; ${params.diamond_resistance_mincov}/100" | bc -l ) ;
  awk -F '\t' '{ if (\$2 != "") { print } }' AMRFinder_complete.tsv > AMRFinder_resistance-only.tsv ;
  """
}
