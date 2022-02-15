process DIGIS {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "${prefix}.gff") null
    else "$filename"
  }
  tag "${prefix}"
  label 'misc'

  input:
  tuple val(prefix), file(genome), file(genbank)

  output:
  // Grab results
  file("digIS")
  tuple val(prefix), file("digIS/results/*.gff")
  tuple val(prefix), file("${prefix}.gff"), file("digIS/results/fastas/${prefix}_is.fa"), file("digIS/results/fastas/${prefix}_is.faa")

  script:
  """
  # activate env
  source activate digIS

  # run digIS
  python3 \$(which digIS_search.py) -i $genome -g $genbank -o digIS

  # deactivate env
  conda deactivate

  # parse digIS to get nucleotide and aminoacide
  # also put ids in uppercase
  # required for annotation merging and sqldb

  ## dir for fastas
  mkdir -p digIS/results/fastas ;

  ## save info in gff
  sed \
    -e 's/id=/ID=/g' \
    digIS/results/${prefix}.gff > ${prefix}.gff ;

  ## get nucl sequences
  gff-toolbox \
    convert \
    -i ${prefix}.gff \
    -f fasta-nt \
    --fasta $genome \
    --fasta_features transposable_element > digIS/results/fastas/${prefix}_is.fa  ;
  
  ## get prot sequences
  gff-toolbox \
    convert \
    -i ${prefix}.gff \
    -f fasta-aa \
    --fasta $genome \
    --fasta_features transposable_element > digIS/results/fastas/${prefix}_is.faa ;
  """
}
