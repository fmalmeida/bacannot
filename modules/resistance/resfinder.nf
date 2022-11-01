process RESFINDER {
    publishDir "${params.output}/${prefix}/resistance", mode: 'copy'
    tag "${prefix}"
    label = [ 'misc', 'process_medium' ]

    input:
    tuple val(prefix), file(genome), val(resfinder_species)
    file(bacannot_db)

    output:
    tuple val(prefix), path("resfinder/ResFinder_results_tab.txt"), emit: results
    tuple val(prefix), path("resfinder/PointFinder_results.txt")  , emit: pointfinder_results
    tuple val(prefix), path("resfinder/args_pheno_table.txt")     , emit: pheno_table
    tuple val(prefix), path("resfinder/results_tab.gff")          , emit: gff
    path("resfinder/*")                                           , emit: all

    when:
    (resfinder_species && resfinder_species != "missing_resfinder")

    script:
    resistance_minid  = params.blast_resistance_minid / 100.00
    resistance_mincov = params.blast_resistance_mincov / 100.00
    if (resfinder_species.toLowerCase() != "other")
    """
    # activate env
    source activate resfinder

    # Make databases available
    ln -rs ${bacannot_db}/resfinder_db/db_* \$(dirname \$(which run_resfinder.py))

    # Run resfinder acquired resistance
    run_resfinder.py \\
        --inputfasta $genome \\
        -o resfinder \\
        --species \"${resfinder_species}\" \\
        --min_cov  ${resistance_mincov} \\
        --threshold ${resistance_minid} \\
        --acquired ;

    # Fix name of pheno table
    mv resfinder/pheno_table.txt resfinder/args_pheno_table.txt &> /dev/null ;

    # Run resfinder pointfinder resistance
    run_resfinder.py \\
        --inputfasta $genome \\
        -o resfinder \\
        --species \"${resfinder_species}\" \\
        --min_cov  ${resistance_mincov} \\
        --threshold ${resistance_minid} \\
        --point ;

    # Fix name of pheno table
    mv resfinder/pheno_table.txt resfinder/mutation_pheno_table.txt &> /dev/null ;

    # Convert to GFF
    resfinder2gff.py \\
        -i resfinder/ResFinder_results_tab.txt > resfinder/results_tab.gff ;
    """

    else if (resfinder_species.toLowerCase() == "other")
    """
    # activate env
    source activate resfinder
    
    # Make databases available
    ln -rs ${bacannot_db}/resfinder_db/db_* \$(dirname \$(which run_resfinder.py))

    # Run resfinder acquired resistance
    run_resfinder.py \\
        --inputfasta $genome \\
        -o resfinder \\
        --species \"${resfinder_species}\" \\
        --min_cov  ${resistance_mincov} \\
        --threshold ${resistance_minid} \\
        --acquired ;

    # Fix name of pheno table
    mv resfinder/pheno_table.txt resfinder/args_pheno_table.txt &> /dev/null ;

    # touch pointfinder
    touch resfinder/PointFinder_results.txt ;

    # Convert to GFF
    resfinder2gff.py \\
        -i resfinder/ResFinder_results_tab.txt > resfinder/results_tab.gff ;
    """
}
