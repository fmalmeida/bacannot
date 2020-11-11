process report {
  publishDir "${params.outdir}/${prefix}/report_files", mode: 'copy'
  label 'renv'
  tag "Rendering HTML reports for virulence, ICEs and AMR genes"

  input:
  tuple val(prefix), file(gff), file(draft), file("prokka_gff"), file(mlst), file(barrnap),
        file(gc_bedGraph), file(gc_chrSizes), file(kofamscan), file(vfdb_blastn), file(victors_blastp),
        file(amrfinder), file(rgi), file(iceberg_blastp), file(phast_blastp), file(phigaro_txt),
        file(genomic_islands), file("methylation"), file("chr.sizes"),
        file(rgi_perfect), file(rgi_strict), file(argminer_out), file(iceberg_blastn), file(plasmids_tsv),
        file(resfinder_tab), file(resfinder_point)

  output:
  file '*.html'

  script:
  """
  cp /work/rscripts/reports/* . ;

  ## Generate Resistance Report
  Rscript -e 'rmarkdown::render("report_resistance.Rmd", params = list(\
    blast_id = ${params.blast_resistance_minid} , \
    blast_cov = ${params.blast_resistance_mincov}, \
    amrfinder = "$amrfinder", \
    query = "${prefix}", \
    rgitool = "$rgi", \
    rgiperfect = "$rgi_perfect", \
    rgistrict = "$rgi_strict", \
    argminer_blastp = "$argminer_out", \
    resfinder_tab = "$resfinder_tab", \
    resfinder_pointfinder = "$resfinder_point", \
    gff = "$gff"))'

  ## Generate Virulence Report
  Rscript -e 'rmarkdown::render("report_virulence.Rmd" , \
  params = list( blast_id = ${params.blast_virulence_minid} , \
                 blast_cov = ${params.blast_virulence_mincov}, \
                 vfdb_blast = "$vfdb_blastn", \
                 gff = "$gff", \
                 victors_blast = "$victors_blastp", \
                 query = "${prefix}"))'

  ## Generate MGEs report
  Rscript -e 'rmarkdown::render("report_MGEs.Rmd", \
  params = list( blast_id = ${params.blast_MGEs_minid}, \
                 blast_cov = ${params.blast_MGEs_mincov}, \
                 phigaro_dir = "${params.outdir}/prophages/phigaro", \
                 phigaro_txt = "$phigaro_txt", \
                 ice_prot_blast = "$iceberg_blastp", \
                 ice_genome_blast = "$iceberg_blastn", \
                 plasmid_finder_tab = "$plasmids_tsv", \
                 query = "${prefix}", \
                 gff = "$gff", \
                 phast_prot_blast = "$phast_blastp" ))'
  """
}
