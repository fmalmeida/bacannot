process RESFINDER {
  publishDir "${params.output}/${prefix}/resistance", mode: 'copy'
  tag "${prefix}"
  label = [ 'misc', 'process_medium' ]

  input:
  tuple val(prefix), file(genome), val(resfinder_species)
  file(bacannot_db)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("resfinder/ResFinder_results_tab.txt")
  tuple val(prefix), file("resfinder/PointFinder_results.txt")
  tuple val(prefix), file("resfinder/args_pheno_table.txt")
  tuple val(prefix), file("resfinder/results_tab.gff")
  file("resfinder/*") // Grab everything

  when:
  (resfinder_species && resfinder_species != "missing_resfinder")

  script:
  resistance_minid  = params.blast_resistance_minid / 100.00
  resistance_mincov = params.blast_resistance_mincov / 100.00
  if (resfinder_species.toLowerCase() != "other")
  """
  # Run resfinder acquired resistance
  run_resfinder.py \\
      --inputfasta $genome \\
      -o resfinder \\
      --species \"${resfinder_species}\" \\
      --min_cov  ${resistance_mincov} \\
      --threshold ${resistance_minid} \\
      --db_path_res ${bacannot_db}/resfinder_db/db_resfinder --acquired || true ;

  # Fix name of pheno table
  mv resfinder/pheno_table.txt resfinder/args_pheno_table.txt &> /dev/null || true ;

  # Run resfinder pointfinder resistance
  run_resfinder.py \\
      --inputfasta $genome \\
      -o resfinder \\
      --species \"${resfinder_species}\" \\
      --min_cov  ${resistance_mincov} \\
      --threshold ${resistance_minid} \\
      --db_path_point ${bacannot_db}/resfinder_db/db_pointfinder --point || true ;

  # Fix name of pheno table
  mv resfinder/pheno_table.txt resfinder/mutation_pheno_table.txt &> /dev/null || true ;

  # Convert to GFF
  resfinder2gff.py \\
      -i resfinder/ResFinder_results_tab.txt > resfinder/results_tab.gff ;
  """

  else if (resfinder_species.toLowerCase() == "other")
  """
  # Run resfinder acquired resistance
  run_resfinder.py \\
      --inputfasta $genome \\
      -o resfinder \\
      --species \"${resfinder_species}\" \\
      --min_cov  ${resistance_mincov} \\
      --threshold ${resistance_minid} \\
      --db_path_res ${bacannot_db}/resfinder_db/db_resfinder --acquired || true ;

  # Fix name of pheno table
  mv resfinder/pheno_table.txt resfinder/args_pheno_table.txt &> /dev/null || true ;

  # touch pointfinder
  touch resfinder/PointFinder_results.txt ;

  # Convert to GFF
  resfinder2gff.py \\
      -i resfinder/ResFinder_results_tab.txt > resfinder/results_tab.gff ;
  """
}
