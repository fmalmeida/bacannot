process resfinder {
  publishDir "${params.outdir}/${prefix}/resistance", mode: 'copy'
  tag "Scanning AMR genes with resfinder"
  label 'main'

  input:
  tuple val(prefix), file(genome), val(resfinder_species)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("resfinder/ResFinder_results_tab.txt")
  tuple val(prefix), file("resfinder/PointFinder_results.txt")
  tuple val(prefix), file("resfinder/args_pheno_table.txt")
  tuple val(prefix), file("resfinder/ResFinder_results_tab.gff")
  file("resfinder/*") // Grab everything

  when:
  (resfinder_species != "missing_resfinder")

  script:
  (resfinder_species.toLowerCase() != "other")
  """
  # Run resfinder acquired resistance
  /work/resfinder/run_resfinder.py --inputfasta $genome -o resfinder --species \"${resfinder_species}\" \
  --min_cov \$(echo "scale=2; ${params.blast_resistance_mincov}/100" | bc -l ) \
  --threshold \$(echo "scale=2; ${params.blast_resistance_minid}/100" | bc -l ) \
  --db_path_res /work/resfinder/db_resfinder --acquired ;

  # Fix name of pheno table
  mv resfinder/pheno_table.txt resfinder/args_pheno_table.txt ;

  # Run resfinder pointfinder resistance
  # this does not run with panel "other", thus we have to create blank files when it fails
  /work/resfinder/run_resfinder.py --inputfasta $genome -o resfinder --species \"${resfinder_species}\" \
  --min_cov \$(echo "scale=2; ${params.blast_resistance_mincov}/100" | bc -l ) \
  --threshold \$(echo "scale=2; ${params.blast_resistance_minid}/100" | bc -l ) \
  --db_path_point /work/resfinder/db_pointfinder --point || touch resfinder/PointFinder_results.txt resfinder/pheno_table.txt ;

  # Fix name of pheno table
  mv resfinder/pheno_table.txt resfinder/mutation_pheno_table.txt &> /dev/null || true ;

  # Convert to GFF
  resfinder2gff.py -i resfinder/ResFinder_results_tab.txt > resfinder/ResFinder_results_tab.gff ;
  """

  (resfinder_species.toLowerCase() == "other")
  """
  # Run resfinder acquired resistance
  /work/resfinder/run_resfinder.py --inputfasta $genome -o resfinder \
  --min_cov \$(echo "scale=2; ${params.blast_resistance_mincov}/100" | bc -l ) \
  --threshold \$(echo "scale=2; ${params.blast_resistance_minid}/100" | bc -l ) \
  --db_path_res /work/resfinder/db_resfinder --acquired || true ;

  # Fix name of pheno table
  mv resfinder/pheno_table.txt resfinder/args_pheno_table.txt ;

  # touch pointfinder
  touch resfinder/PointFinder_results.txt ;

  # Convert to GFF
  resfinder2gff.py -i resfinder/results_tab.txt > resfinder/results_tab.gff ;
  """
}
