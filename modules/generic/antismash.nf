process antismash {
  publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "$filename"
  }
  tag "Scanning genome with antiSMASH"
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
  antismash --output-dir antiSMASH --genefinding-tool none -c ${params.threads} $genbank

  # convert results to gff
  for gbk in antiSMASH/*.region*.gbk ; do
    seqret -sequence \${gbk} -feature -fformat genbank -fopenfile \${gbk} -osformat gff -osname_outseq \${gbk%%.gbk} -auto ;
  done
  """
}
