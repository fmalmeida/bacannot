process ANTISMASH {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "$filename"
  }
  tag "${prefix}"
  label = [ 'misc', 'process_medium' ]

  input:
  tuple val(prefix), file(genbank)
  file(bacannot_db)

  output:
  // Grab results
  tuple val(prefix), path("antiSMASH/regions.gff")
  path("antiSMASH")
  path("*_version.txt")

  script:
  """
  # Activate env
  export PATH=/opt/conda/envs/antismash/bin:\$PATH
  
  # Get tool version
  antismash --version > antismash_version.txt ;

  # Run tool
  antismash \\
    --output-dir antiSMASH \\
    --genefinding-tool none \\
    -c ${params.threads} \\
    --databases ${bacannot_db}/antismash_db \\
    $genbank ;

  # enter results dir
  cd antiSMASH ;

  # produce gff from main results
  genbank="${genbank}"
  seqret \\
    -sequence \${genbank} \\
    -feature \\
    -fformat genbank \\
    -fopenfile \${genbank} \\
    -osformat gff \\
    -osname_outseq \${genbank%%.gbk} \\
    -auto ;

  # get the locus tags annotated as list
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
    \${genbank%%.gbk}.gff > regions.gff ;
  """
}
