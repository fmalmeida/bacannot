process kofamscan {
  publishDir "${params.outdir}", mode: 'copy'
  container = 'fmalmeida/bacannot:kofamscan'
  tag "Executing KOfamscan - Its output can be directly put in KEGG Mapper for visualization"

  input:
  file 'proteins.faa'

  output:
  file "KOfamscan/${params.prefix}_ko_*" // Get all files to input directory
  file "KOfamscan/${params.prefix}_ko_forKEGGMapper.txt" // Kegg-mapper file

  script:
  """
  mkdir KOfamscan ;
  kofamscan -o KOfamscan/${params.prefix}_ko_detailed.txt --cpu=${params.threads} proteins.faa ;
  kofamscan -o KOfamscan/${params.prefix}_ko_forKEGGMapper.txt --cpu=${params.threads} -f mapper-one-line proteins.faa ;
  """
}
