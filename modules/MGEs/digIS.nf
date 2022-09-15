process DIGIS {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "${prefix}.gff") null
    else if (filename == "${prefix}_IS.gff") null
    else "$filename"
  }
  tag "${prefix}"
  label = [ 'misc', 'process_low' ]

  input:
  tuple val(prefix), path(genome), path(genbank)

  output:
  path("digIS")                                         , emit: all
  tuple val(prefix), path("digIS/results/${prefix}.gff"), emit: gff
  tuple val(prefix), path("${prefix}_IS.gff"), path("digIS/results/fastas/${prefix}_IS.fa"), path("digIS/results/fastas/${prefix}_IS.faa"), emit: gff_and_sequences

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
  sed \\
    -e 's/id=/ID=/g' \\
    digIS/results/${prefix}.gff > ${prefix}_IS.gff ;

  ## get nucl sequences
  gff-toolbox \\
    convert \\
    -i ${prefix}_IS.gff  \\
    -f fasta-nt \\
    --fasta $genome \\
    --fasta_features transposable_element > digIS/results/fastas/${prefix}_IS.fa  ;
  
  ## get prot sequences
  gff-toolbox \\
    convert \\
    -i ${prefix}_IS.gff  \\
    -f fasta-aa \\
    --fasta $genome \\
    --fasta_features transposable_element > digIS/results/fastas/${prefix}_IS.faa ;
  """
}
