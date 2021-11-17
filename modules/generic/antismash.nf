process antismash {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "$filename"
  }
  tag "${prefix}"
  label 'smash'

  input:
  tuple val(prefix), file(genbank)

  output:
  // Grab results
  tuple val(prefix), file("antiSMASH")
  file("*_version.txt")

  script:
  """
  # activate env
  source activate antismash ;

  # Get tool version
  antismash --version > antismash_version.txt ;

  # Run tool
  antismash --output-dir antiSMASH --genefinding-tool none -c ${params.threads} $genbank ;

  # enter results dir
  cd antiSMASH ;

  # produce gff from main results
  genbank="${genbank}"
  seqret -sequence \${genbank} -feature -fformat genbank -fopenfile \${genbank} -osformat gff -osname_outseq \${genbank%%.gbk} -auto ;

  # get the locus tags annotated as list
  grep "locus_tag" *region*gbk | cut -f 2 -d "=" | tr -d '"' | sort -u > gene_ids.lst ;

  # subset regions GFF from main GFF for JBrowse
  grep -w -f gene_ids.lst \${genbank%%.gbk}.gff > regions.gff ;
  """
}
