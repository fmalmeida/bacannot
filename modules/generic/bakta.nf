process BAKTA {
    publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
        if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
        else if (filename == "annotation") "$filename"
        else null
    }
    tag "${prefix}"
    label = [ 'process_medium', 'error_retry' ]

    conda "bioconda::bakta=1.7.0"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bakta:1.7.0--pyhdfd78af_1' :
        'quay.io/biocontainers/bakta:1.7.0--pyhdfd78af_1' }"

    input:
    tuple val(prefix), val(entrypoint), file(sread1), file(sread2), file(sreads), file(lreads), val(lr_type), file(fast5), file(assembly), val(resfinder_species)
    file(bakta_db)

    output:
    // Grab all outputs
    tuple val(prefix), path("annotation"), emit: all
    // Outputs must be linked to each prefix (tag)
    tuple val(prefix), path("annotation/${prefix}.gff3"), emit: gff // annotation in gff format
    tuple val(prefix), path("annotation/${prefix}.gbff"), emit: gbk // annotation in gbk format
    tuple val(prefix), path("annotation/${prefix}.fna") , emit: genome // renamed genome
    tuple val(prefix), path("annotation/${prefix}.faa") , emit: proteins // gene aa sequences
    tuple val(prefix), path("annotation/${prefix}.ffn") , emit: genes // gene nt sequences
    tuple val(prefix), path("annotation/${prefix}.fna"), path("${lreads}"), path("${fast5}"), emit: genome_with_fast5 // For methylation calling
    tuple val(prefix), path("annotation/${prefix}.fna"), val("${resfinder_species}"), emit: genome_with_species // For resfinder
    tuple val(prefix), path("annotation/${prefix}.txt") , emit: summary // bakta stats
    path('bakta_version.txt'), emit: version // Save bakta version

    script:
    """
    # download amrfinder db if not available
    [ ! -d ${bakta_db}/amrfinderplus-db ] && amrfinder_update --database ${bakta_db}/amrfinderplus-db

    # Save bakta version
    bakta --version &> bakta_version.txt ;

    # clean headers char limit
    awk '{ if (\$0 ~ />/) print substr(\$0,1,21) ; else print \$0 }' $assembly > cleaned_header.fasta

    # Run bakta
    bakta \\
        --output annotation \\
        --threads $task.cpus \\
        --min-contig-length 200 \\
        --prefix ${prefix} \\
        --strain '${prefix}' \\
        --db $bakta_db \\
        cleaned_header.fasta
    
    # fix fasta headers
    cut -f 1 -d ' ' annotation/${prefix}.fna > tmp.fa
    cat tmp.fa > annotation/${prefix}.fna
    rm tmp.fa
    """
}
