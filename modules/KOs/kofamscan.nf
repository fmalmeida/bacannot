process KOFAMSCAN {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "$filename"
  }
  errorStrategy 'retry'
  maxRetries 2
  tag "${prefix}"
  label 'kofam'

  input:
  tuple val(prefix), file('proteins.faa')

  output:
  // Grab all outputs
  file("KOfamscan") // Get all files to input directory
  tuple val(prefix), file("KOfamscan/${prefix}_ko_forKEGGMapper.txt") // Kegg-mapper file

  script:
  """
  # Get kofamscan version
  kofamscan -v > kofamscan_version.txt

  # Create dir for results
  mkdir KOfamscan ;

  # Run kofamscan with detailed output
  kofamscan -o KOfamscan/${prefix}_ko_detailed.txt --keep-tabular --cpu=${params.threads} proteins.faa ;

  # Re-run kofamscan with mapper-output
  kofamscan -o KOfamscan/${prefix}_ko_forKEGGMapper.txt --reannotate --cpu=${params.threads} -f mapper proteins.faa ;
  """
}
