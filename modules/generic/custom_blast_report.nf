process custom_blast_report {
  publishDir "${params.output}/${prefix}/report_files/custom_databases", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf(".html") > 0) "report_${customDB}.html"
    else "$filename"
  }
  label 'renv'
  tag "Rendering HTML reports for the custom db annotations"

  input:
  tuple val(prefix), val(customDB), file(custom_blast), file(custom_gff)

  output:
  file '*.html'

  script:
  """
  #!/usr/bin/env Rscript

  ## Copy reports
  system("cp /work/reports/* .") ;

  ## Remove empty files
  system("rm -f input.??") ;
  system("rm -f input.?") ;

  ## Generate Resistance Report
  rmarkdown::render("report_custom_blast.Rmd", params = list(\
    blast_id = ${params.blast_custom_minid} , \
    blast_cov = ${params.blast_custom_mincov}, \
    query = "${prefix}", \
    custom_blast = "$custom_blast", \
    blast_db = "${customDB}", \
    blast_gff = "$custom_gff")) ;
  """
}
