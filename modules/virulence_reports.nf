process virulence_report {
  publishDir "${params.outdir}/${prefix}/report_files", mode: 'copy'
  label 'reports'

  input:
  tuple val(prefix), file(gff), file(draft), file("prokka_gff"), file(mlst), file(barrnap),
        file(gc_bedGraph), file(gc_chrSizes), file(kofamscan), file(vfdb), file(victors),
        file(amrfinder), file(rgi), file(iceberg), file(phast), file(phigaro),
        file(genomic_islands), file("cpg"), file("gpc"), file("dam"), file("dcm"),
        file("chr.sizes"), file(phigaro_txt)

  output:
  file '*.html'

  script:
  """
  cp /work/rscripts/reports/* . ;

  ## Generate Virulence Report
  Rscript -e 'rmarkdown::render("report_virulence.Rmd" , \
  params = list( vfdb_blast = "$vfdb", \
                 blast_id = ${params.diamond_virulence_identity} , \
                 blast_cov = ${params.diamond_virulence_queryCoverage}, \
                 gff = "$gff", \
                 victors_blast = "$victors", \
                 query = "${prefix}"))'
  """
}
