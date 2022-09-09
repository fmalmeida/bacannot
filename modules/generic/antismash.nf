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
  def gbk_suffix = (params.bakta_db) ? "gbff" : "gbk"
  def gbk_prefix = "${genbank.baseName}" - "${gbk_suffix}"
  """
  # Activate env
  export PATH=/opt/conda/envs/antismash/bin:\$PATH
  
  # Get tool version
  antismash --version > antismash_version.txt ;

  # Run tool
  antismash \\
    --output-dir antiSMASH \\
    --genefinding-tool none \\
    -c $task.cpus \\
    --databases ${bacannot_db}/antismash_db \\
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
  """
}
