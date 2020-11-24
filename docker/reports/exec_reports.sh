## Generate Resistance Report
Rscript -e 'rmarkdown::render("report_resistance.Rmd", params = list(
  blast_id = 90 ,
  blast_cov = 90,
  amrfinder = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/resistance/AMRFinderPlus/AMRFinder_resistance-only.tsv",
  query = "KpV3",
  rgitool = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/resistance/RGI/RGI_KpV3.txt",
  rgiperfect = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/resistance/RGI/Perfect_RGI_KpV3_hits.txt",
  rgistrict = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/resistance/RGI/Strict_RGI_KpV3_hits.txt",
  rgi_heatmap = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/resistance/RGI/heatmap/RGI_heatmap-1.png",
  argminer_blastp = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/resistance/ARGMiner/KpV3_argminer_blastp_onGenes.summary.txt",
  resfinder_tab = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/resistance/resfinder/results_tab.txt",
  resfinder_pointfinder = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/resistance/resfinder/PointFinder_results.txt",
  resfinder_phenotype = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/resistance/resfinder/pheno_table.txt",
  gff = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/gffs/KpV3.gff"))'

## Generate Virulence Report
Rscript -e 'rmarkdown::render("report_virulence.Rmd" , 
  params = list( blast_id = 90, 
                 blast_cov = 90, 
                 vfdb_blast = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/virulence/vfdb/KpV3_vfdb_blastn_onGenes.txt", 
                 gff = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/gffs/KpV3.gff", 
                 victors_blast = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/virulence/victors/KpV3_victors_blastp_onGenes.txt", 
                 query = "KpV3"))'
                 
## Generate MGEs report
Rscript -e 'rmarkdown::render("report_MGEs.Rmd", 
  params = list( blast_id = 65, 
                 blast_cov = 65, 
                 phigaro_dir = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/prophages/phigaro", 
                 phigaro_txt = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/prophages/phigaro/KpV3_phigaro.tsv",
                 phispy_tsv = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/prophages/PhiSpy/prophage.tsv",
                 ice_prot_blast = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/ICEs/KpV3_iceberg_blastp_onGenes.txt", 
                 ice_genome_blast = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/ICEs/KpV3_iceberg_blastn_onGenome.summary.txt", 
                 plasmid_finder_tab = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/plasmids/plasmidfinder/results_tab.tsv",
                 platon_tsv = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/plasmids/platon/KpV3.tsv",
                 query = "KpV3",
                 gi_image = "/Volumes/falmeida1TB/bacannot_teste/output/dataset1/genomic_islands/plots/product_label/contig_1_GI_1.png",
                 gff = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/gffs/KpV3.gff", 
                 phast_prot_blast = "/Volumes/falmeida1TB/bacannot_teste/output/KpV3/prophages/phast_db/KpV3_phast_blastp_onGenes.txt" ))'
                 
## Clean dir
rm -f ._*
                 
                 