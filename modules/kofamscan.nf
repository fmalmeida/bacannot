process kofamscan {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  tag "Executing KOfamscan - Its outputs can be viewed in KEGG-mapper"
  label 'kofam'

  input:
  tuple val(prefix), file('proteins.faa')

  output:
  // Grab all outputs
  file("KOfamscan/*") // Get all files to input directory
  tuple val(prefix), file("KOfamscan/${prefix}_ko_forKEGGMapper.txt") // Kegg-mapper file

  script:
  """
  # Create dir for results
  mkdir KOfamscan ;

  # Run kofamscan with detailed output
  kofamscan -o KOfamscan/${prefix}_ko_detailed.txt --cpu=${params.threads} proteins.faa ;

  # Run kofamscan with mapper-output
  kofamscan -o KOfamscan/${prefix}_ko_forKEGGMapper.txt --cpu=${params.threads} -f mapper proteins.faa ;
  """
}
