process PROKKA {
    publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
        if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
        else if (filename == "annotation") "$filename"
        else null
    }
    tag "${prefix}"
    label = [ 'process_medium' ]

    conda "bioconda::prokka=1.14.6"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl5321hdfd78af_4' :
        'quay.io/biocontainers/prokka:1.14.6--pl5321hdfd78af_4' }"

    input:
    tuple val(prefix), val(entrypoint), file(sread1), file(sread2), file(sreads), file(lreads), val(lr_type), file(fast5), file(assembly), val(resfinder_species)
    file(bacannot_db)

    output:
    // Grab all outputs
    tuple val(prefix), path("annotation"), emit: all
    // Outputs must be linked to each prefix (tag)
    tuple val(prefix), path("annotation/${prefix}.gff"), emit: gff
    tuple val(prefix), path("annotation/${prefix}.gbk"), emit: gbk
    tuple val(prefix), path("annotation/${prefix}.fna"), emit: genome
    tuple val(prefix), path("annotation/${prefix}.faa"), emit: proteins
    tuple val(prefix), path("annotation/${prefix}.ffn"), emit: genes
    tuple val(prefix), path("annotation/${prefix}.fna"), path("${lreads}"), path("${fast5}"), emit: genome_with_fast5
    tuple val(prefix), path("annotation/${prefix}.fna"), val("${resfinder_species}"), emit: genome_with_species
    tuple val(prefix), path("annotation/${prefix}.txt"), emit: summary
    path('prokka_version.txt'), emit: version

    script:
    kingdom = (params.prokka_kingdom)      ? "--kingdom ${params.prokka_kingdom}"    : ''
    gcode   = (params.prokka_genetic_code) ? "--gcode ${params.prokka_genetic_code}" : ''
    rnammer = (params.prokka_use_rnammer)  ? "--rnammer"                             : ''
    models  = (params.prokka_use_pgap)     ? "PGAP_NCBI.hmm"                         : "TIGRFAMs_15.0.hmm"
    """
    #!/usr/bin/env bash
    
    # save prokka version
    prokka -v &> prokka_version.txt ;

    # where are default prokka dbs?
    dbs_dir=\$(prokka --listdb 2>&1 >/dev/null |  grep "databases in" | cut -f 4 -d ":" | tr -d " ") ;

    # get hmms that shall be used
    # PGAP contains TIGRFAM hmm models. When not skipping PGAP, TIGRFAM is not loaded.
    cp -r \$dbs_dir prokka_db
    cp ${bacannot_db}/prokka_db/${models} prokka_db/hmm

    # hmmpress
    ( cd  prokka_db/hmm/ ; for i in *.hmm ; do hmmpress -f \$i ; done )

    # clean headers char limit
    awk '{ if (\$0 ~ />/) print substr(\$0,1,21) ; else print \$0 }' $assembly > cleaned_header.fasta

    # run prokka
    prokka \\
        --dbdir prokka_db \\
        $kingdom $gcode $rnammer \\
        --outdir annotation \\
        --cpus $task.cpus \\
        --mincontiglen 200 \\
        --prefix ${prefix} \\
        --genus '' \\
        --species '' \\
        --strain \"${prefix}\" \\
        cleaned_header.fasta
    
    # remove tmp dir to gain space
    rm -r prokka_db
    """
}
