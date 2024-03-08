process ANTISMASH {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "$filename"
  }
  tag "${prefix}"
  label = [ 'misc', 'process_medium' ]

  // if (params.running_engine = 'singularity') { runOptions = '--writable-tmpfs -e --no-home -B $PWD' }

  input:
  tuple val(prefix), file(genbank)
  file(bacannot_db)

  output:
  tuple val(prefix), path("antiSMASH/regions.gff"), emit: gff, optional: true
  tuple val(prefix), path("antiSMASH")            , emit: all
  path("*_version.txt")                           , emit: version

  script:
  def gbk_suffix = (params.bakta_db) ? "gbff" : "gbk"
  def gbk_prefix = "${genbank.baseName}" - "${gbk_suffix}"
  def antismash_version='6.1.1'

  if (params.running_engine == 'singularity')
  """  
  # Get tool version
  antismash --version > antismash_version.txt ;

  # activate env
  mkdir local-install
  export PYTHONUSERBASE=./local-install
  export PATH=/opt/conda/envs/antismash/bin:\$PATH

  # singularity has many read-write permissions for this tool
  wget https://dl.secondarymetabolites.org/releases/${antismash_version}/antismash-${antismash_version}.tar.gz
  tar zxvf antismash-${antismash_version}.tar.gz
  python -m pip install --user ./antismash-${antismash_version}
  export PYTHONPATH=\$(realpath \$( find ./local-install -name 'site-packages' ))

  # Run tool
  ./local-install/bin/antismash \\
    --output-dir antiSMASH \\
    --genefinding-tool none \\
    --databases ${bacannot_db}/antismash_db \\
    -c $task.cpus \\
    $genbank ;
  
  # enter results dir
  cd antiSMASH ;

  # produce gff from main results
  seqret \\
    -sequence ${gbk_prefix}.gbk \\
    -feature \\
    -fformat genbank \\
    -fopenfile ${gbk_prefix}.gbk \\
    -osformat gff \\
    -osname_outseq ${gbk_prefix} \\
    -auto ;

  # get the locus tags annotated as list
  # only when results exist
  if ls *region*gbk 1> /dev/null 2>&1; then

    grep \\
      "locus_tag" \\
      *region*gbk | \\
      cut \\
      -f 2 \\
      -d "=" | \\
      tr -d '"' | \\
      sort -u > gene_ids.lst ;

    # subset regions GFF from main GFF for JBrowse
    grep \\
      -w \\
      -f gene_ids.lst \\
      ${gbk_prefix}.gff > regions.gff ;
  
  fi
  """

  else
  """  
  # Get tool version
  antismash --version > antismash_version.txt ;

  # Run tool
  antismash \\
    --output-dir antiSMASH \\
    --genefinding-tool none \\
    --databases ${bacannot_db}/antismash_db \\
    -c $task.cpus \\
    $genbank ;
  
  # enter results dir
  cd antiSMASH ;

  # produce gff from main results
  seqret \\
    -sequence ${gbk_prefix}.gbk \\
    -feature \\
    -fformat genbank \\
    -fopenfile ${gbk_prefix}.gbk \\
    -osformat gff \\
    -osname_outseq ${gbk_prefix} \\
    -auto ;

  # get the locus tags annotated as list
  # only when results exist
  if ls *region*gbk 1> /dev/null 2>&1; then

    grep \\
      "locus_tag" \\
      *region*gbk | \\
      cut \\
      -f 2 \\
      -d "=" | \\
      tr -d '"' | \\
      sort -u > gene_ids.lst ;

    # subset regions GFF from main GFF for JBrowse
    grep \\
      -w \\
      -f gene_ids.lst \\
      ${gbk_prefix}.gff > regions.gff ;
  
  fi
  """

}
