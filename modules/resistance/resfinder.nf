process resfinder {
  publishDir "${params.outdir}/${prefix}/resistance", mode: 'copy'
  tag "Scanning AMR genes with resfinder"
  label 'main'

  input:
  tuple val(prefix), file(genome), val(resfinder_species)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("resfinder/results_tab.txt")
  tuple val(prefix), file("resfinder/PointFinder_results.txt")
  tuple val(prefix), file("resfinder/args_pheno_table.txt")
  file("resfinder") // Grab everything

  when:
  (resfinder_species != "missing_resfinder")

  script:
  """
  # Activate env
  source activate Resfinder;

  # Run resfinder acquired resistance
  /work/resfinder/run_resfinder.py --inputfasta $genome -o resfinder --species \"${resfinder_species}\" \
  --min_cov \$(echo "scale=2; ${params.blast_resistance_mincov}/100" | bc -l ) \
  --threshold \$(echo "scale=2; ${params.blast_resistance_minid}/100" | bc -l ) \
  --db_path_res /work/resfinder/db_resfinder --acquired || true ;

  # Fix name of pheno table
  mv resfinder/pheno_table.txt resfinder/args_pheno_table.txt

  # Run resfinder pointfinder resistance
  /work/resfinder/run_resfinder.py --inputfasta $genome -o resfinder --species \"${resfinder_species}\" \
  --min_cov \$(echo "scale=2; ${params.blast_resistance_mincov}/100" | bc -l ) \
  --threshold \$(echo "scale=2; ${params.blast_resistance_minid}/100" | bc -l ) \
  --db_path_point /work/resfinder/db_pointfinder --point || true ;

  # Fix name of pheno table
  mv resfinder/pheno_table.txt resfinder/mutation_pheno_table.txt &> /dev/null || true ;
  """
}
