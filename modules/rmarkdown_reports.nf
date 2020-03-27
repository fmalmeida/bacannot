process report {
  publishDir "${outDir}/report_files", mode: 'copy'
  container 'fmalmeida/bacannot:renv'

  input:
  val x from finish
  file 'final.gff' from final_gff
  file rgi_table from rgi_output
  file rgi_perfect from rgi_perfect
  file rgi_strict from rgi_strict
  file amrfinder_result from amrfinder_output
  file amrfinder_summary from amrfinderplus_table
  file 'resistance.gff' from resistance_gff
  file vfdb_blast from vfdb_blast_genes2
  file victors_blast from victors_blast_genes2
  file ice_blast from ice_blast_genes2
  //file virsorter_csv from virsorter_csv.ifEmpty('virsorter_empty')
  file phigaro_txt from phigaro_txt.ifEmpty('phigaro_empty')
  file phast_blast from phast_blast_genes2

  output:
  file '*.html'

  script:
  """
  cp /work/rscripts/*.Rmd . ;

  ## Generate Resistance Report
  Rscript -e 'rmarkdown::render("report_resistance.Rmd", params = list(\
    amrfinder = "$amrfinder_result", \
    query = "${params.prefix}", \
    rgitool = "$rgi_table", \
    rgiperfect = "$rgi_perfect", \
    rgistrict = "$rgi_strict", \
    gff_resistance = "resistance.gff", \
    ncbi_amr = "${amrfinder_summary}"))'

  ## Generate Virulence Report
  Rscript -e 'rmarkdown::render("report_virulence.Rmd" , \
  params = list( vfdb_blast = "${vfdb_blast}", \
                 blast_id = ${params.diamond_virulence_identity} , \
                 blast_cov = ${params.diamond_virulence_queryCoverage},
                 gff = "final.gff",
                 victors_blast = "${victors_blast}",
                 query = "${params.prefix}"))'

  ## Generate MGEs report
  Rscript -e 'rmarkdown::render("report_MGEs.Rmd", params = list( \
                 phigaro_dir = "../prophages/phigaro",
                 phigaro_txt = "${phigaro_txt}",
                 ice_prot_blast = "${ice_blast}",
                 query = "${params.prefix}",
                 gff = "final.gff",
                 blast_id = ${params.diamond_MGEs_identity},
                 blast_cov = ${params.diamond_MGEs_queryCoverage},
                 phast_blast = "${phast_blast}"))'
  """
}
