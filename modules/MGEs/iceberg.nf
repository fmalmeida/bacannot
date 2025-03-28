process ICEBERG {
    publishDir "${params.output}/${prefix}/ICEs", mode: 'copy'
    tag "${prefix}"
    label = [ 'misc', 'process_low' ]

    input:
    tuple val(prefix), file(genes_aa)
    tuple val(prefix), file(genome)
    file(bacannot_db)

    output:
    tuple val(prefix), path("${prefix}_iceberg_blastp_onGenes.summary.txt") , emit: genes_summary
    tuple val(prefix), path("${prefix}_iceberg_blastp_onGenes.txt")         , emit: results
    tuple val(prefix), path("${prefix}_iceberg_blastn_onGenome.summary.txt"), emit: genome_summary
    tuple val(prefix), path('*.txt')                                        , emit: all

    script:
    """
    # ICEberg is a protein and nucleotide dabatase
    # In protein are the genes found inside ICEs
    # In nucleotide are the full-length ICEs

    ## Checking ICE genes
    ## With predicted gene sequences
    run_blasts.py \\
        blastp \\
        --query $genes_aa \\
        --db ${bacannot_db}/iceberg_db/diamond.dmnd \\
        --minid ${params.blast_MGEs_minid} \\
        --mincov ${params.blast_MGEs_mincov} \\
        --threads $task.cpus \\
        --out ${prefix}_iceberg_blastp_onGenes.txt --2way | \\
    sed -e 's/GENE/ICEBERG_ID/g' -e 's/;//g' > ${prefix}_iceberg_blastp_onGenes.summary.txt ;

    ## Checking for full-length ICEs
    ### The blast db was throwing errors
    makeblastdb \\
        -dbtype nucl \\
        -in ${bacannot_db}/iceberg_db/sequences \\
        -out sequences ;
    run_blasts.py \\
        blastn \\
        --query $genome \\
        --db sequences \\
        --minid 0 \\
        --mincov 0 \\
        --threads $task.cpus \\
        --out ${prefix}_iceberg_blastn_onGenome.txt | \\
    sed -e 's/GENE/ICEBERG_ID/g' -e 's/;//g' > ${prefix}_iceberg_blastn_onGenome.summary.txt ;
    """
}
