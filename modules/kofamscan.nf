process kofamscan {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  container = 'fmalmeida/bacannot:kofamscan'
  tag "Executing KOfamscan - Its output can be directly put in KEGG Mapper for visualization"

  input:
  file 'proteins.faa'
  val(prefix)

  output:
  file "KOfamscan/${prefix}_ko_*" // Get all files to input directory
  file "KOfamscan/${prefix}_ko_forKEGGMapper.txt" // Kegg-mapper file

  script:
  """
  mkdir KOfamscan ;
  kofamscan -o KOfamscan/${prefix}_ko_detailed.txt --cpu=${params.threads} proteins.faa ;
  kofamscan -o KOfamscan/${prefix}_ko_forKEGGMapper.txt --cpu=${params.threads} -f mapper-one-line proteins.faa ;
  """
}
