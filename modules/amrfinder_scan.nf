process amrfinder {
  publishDir "${params.outdir}/${prefix}/resistance/AMRFinderPlus", mode: 'copy'
  container = 'fmalmeida/bacannot:latest'
  tag "Resistance genes annotation with AMRFinderPlus"

  input:
  file proteins
  val(prefix)

  output:
  file "AMRFinder_output.tsv"

  script:
  """
  source activate AMRFINDERPLUS ;
  amrfinder -p $proteins --plus -o AMRFinder_output.tsv
  """
}
