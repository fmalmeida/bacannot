process digis {
  publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "${prefix}.gff") null
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
  tuple val(prefix), file("${prefix}.gff"), file("digIS/results/fastas/${prefix}_is.fa"), file("digIS/results/fastas/${prefix}_is.faa") optional true

  script:
  """
  # activate env
  source activate digIS ;

  # run digIS
  python3 /work/digIS/digIS_search.py -i $genome -g $genbank -o digIS

  # parse digIS to get nucleotide and aminoacide
  # also put ids in uppercase
  # required for annotation merging and sqldb
  conda deactivate ;
  mkdir -p digIS/results/fastas ;
  sed -e 's/id=/ID=/g' digIS/results/${prefix}.gff > ${prefix}.gff ;
  gff-toolbox convert -i ${prefix}.gff -f fasta-nt --fasta $genome --fasta_features transposable_element > digIS/results/fastas/${prefix}_is.fa  ;
  gff-toolbox convert -i ${prefix}.gff -f fasta-aa --fasta $genome --fasta_features transposable_element > digIS/results/fastas/${prefix}_is.faa ;
  """
}
