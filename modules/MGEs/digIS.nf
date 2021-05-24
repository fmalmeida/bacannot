process digis {
  publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "$filename"
  }
  tag "Scanning for Insertion Sequences with digIS"
  label 'main'

  input:
  tuple val(prefix), file(genome), file(genbank)

  output:
  // Grab results
  file("digIS")
  tuple val(prefix), file("digIS/results/*.gff") optional true

  script:
  """
  # activate env
  source activate digIS ;

  # run digIS
  python3 /work/digIS/digIS_search.py -i $genome -g $genbank -o digIS
  """
}
