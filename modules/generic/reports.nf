process REPORT {
  publishDir "${params.output}/${prefix}/report_files", mode: 'copy'
  label = [ 'renv', 'process_medium' ]
  tag "${prefix}"

  input:
  tuple val(prefix), path(prokka_stats), path(gff), path(barrnap), path(mlst), path(keggsvg), path(refseq_masher_txt), path(amrfinder), path(rgi), path(rgi_parsed), path(rgi_heatmap), path(argminer_out), path(resfinder_tab), path(resfinder_point), path(resfinder_phenotable), path(vfdb_blastn), path(victors_blastp), path(phigaro_txt), path(phispy_tsv), path(iceberg_blastp), path(iceberg_blastn), path(plasmids_tsv), path(platon_tsv), path(gi_image), path(phast_blastp), path(digIS)
  
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

  ## Generate generic Report
  rmarkdown::render("report_general.Rmd" , \
  params = list( prokka  = "$prokka_stats", \
                 kegg    = "$keggsvg", \
                 barrnap = "$barrnap", \
                 mlst    = "$mlst", \
                 refseq_masher = "$refseq_masher_txt", \
                 query = "${prefix}")) ;

  ## Generate Resistance Report
  rmarkdown::render("report_resistance.Rmd", params = list(\
    blast_id = ${params.blast_resistance_minid} , \
    blast_cov = ${params.blast_resistance_mincov}, \
    amrfinder = "$amrfinder", \
    query = "${prefix}", \
    rgitool = "$rgi", \
    rgiparsed = "$rgi_parsed", \
    rgi_heatmap = "$rgi_heatmap", \
    argminer_blastp = "$argminer_out", \
    resfinder_tab = "$resfinder_tab", \
    resfinder_pointfinder = "$resfinder_point", \
    resfinder_phenotype = "$resfinder_phenotable", \
    gff = "$gff")) ;

  ## Generate Virulence Report
  rmarkdown::render("report_virulence.Rmd" , \
  params = list( blast_id = ${params.blast_virulence_minid} , \
                 blast_cov = ${params.blast_virulence_mincov}, \
                 vfdb_blast = "$vfdb_blastn", \
                 gff = "$gff", \
                 victors_blast = "$victors_blastp", \
                 query = "${prefix}")) ;

  ## Generate MGEs report
  rmarkdown::render("report_MGEs.Rmd", \
  params = list( blast_id = ${params.blast_MGEs_minid}, \
                 blast_cov = ${params.blast_MGEs_mincov}, \
                 phigaro_dir = "${params.output}/prophages/phigaro", \
                 phigaro_txt = "$phigaro_txt", \
                 phispy_tsv = "$phispy_tsv", \
                 ice_prot_blast = "$iceberg_blastp", \
                 ice_genome_blast = "$iceberg_blastn", \
                 plasmid_finder_tab = "$plasmids_tsv", \
                 platon_tsv = "$platon_tsv", \
                 query = "${prefix}", \
                 gi_image = "$gi_image", \
                 digis = "$digIS", \
                 gff = "$gff", \
                 phast_prot_blast = "$phast_blastp" )) ;
  """
}
