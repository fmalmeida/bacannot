process REPORT {
  publishDir "${params.output}/${prefix}/report_files", mode: 'copy'
  label = [ 'renv', 'process_medium' ]
  tag "${prefix}"

  input:
  tuple val(prefix), file('annotation_stats.tsv'), file(gff), file(barrnap), file(mlst), file(keggsvg), file(refseq_masher_txt), file(amrfinder), file(rgi), file(rgi_parsed), file(rgi_heatmap), file(argminer_out), file(resfinder_tab), file(resfinder_point), file(resfinder_phenotable), file(vfdb_blastn), file(victors_blastp), file(phigaro_txt), file(phispy_tsv), file(iceberg_blastp), file(iceberg_blastn), file(plasmids_tsv), file(platon_tsv), file(gi_image), file(phast_blastp), file(digIS), file(integronfinder)
  
  output:
  path '*.html', emit: results

  script:
  def generic_annotator = (params.bakta_db) ? "bakta" : "prokka"
  """
  #!/usr/bin/env Rscript

  ## Copy reports
  system("cp /work/reports/* .") ;

  ## Remove empty files
  system("rm -f input.??") ;
  system("rm -f input.?") ;

  ## Generate generic Report
  rmarkdown::render("report_general.Rmd" , \
    params = list( 
      generic_annotation  = "annotation_stats.tsv", \
      generic_annotator   = "${generic_annotator}", \
      kegg    = "$keggsvg", \
      barrnap = "$barrnap", \
      mlst    = "$mlst", \
      refseq_masher = "$refseq_masher_txt", \
      query = "${prefix}"
    )
  ) ;

  ## Generate Resistance Report
  rmarkdown::render("report_resistance.Rmd", \
    params = list(\
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
      generic_annotator   = "${generic_annotator}", \
      gff = "$gff"
    )
  ) ;

  ## Generate Virulence Report
  rmarkdown::render("report_virulence.Rmd" , \
    params = list( 
      blast_id = ${params.blast_virulence_minid} , \
      blast_cov = ${params.blast_virulence_mincov}, \
      vfdb_blast = "$vfdb_blastn", \
      gff = "$gff", \
      victors_blast = "$victors_blastp", \
      query = "${prefix}"
    )
  ) ;

  ## Generate MGEs report
  rmarkdown::render("report_MGEs.Rmd", \
    params = list( 
      blast_id = ${params.blast_MGEs_minid}, \
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
      integronfinder = "$integronfinder", \
      gff = "$gff", \
      phast_prot_blast = "$phast_blastp"
    )
  ) ;
  """
}
