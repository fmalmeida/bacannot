process amrfinder {
  publishDir "${params.outdir}/${prefix}/resistance/AMRFinderPlus", mode: 'copy'
  container = 'fmalmeida/bacannot:latest'
  tag "Resistance genes annotation with AMRFinderPlus"

  input:
  file proteins
  val(prefix)

  output:
  file "AMRFinder_resistance-only.tsv"
  file "AMRFinder_complete.tsv"

  script:
  """
  source activate AMRFINDERPLUS ;
  amrfinder -p $proteins --plus -o AMRFinder_complete.tsv ;
  awk -F '\t' '{ if (\$2 != "") { print } }' AMRFinder_complete.tsv > AMRFinder_resistance-only.tsv ;
  """
}
