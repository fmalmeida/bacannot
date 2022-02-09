process PROKKA {
    publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
      if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
      else if (filename == "annotation") "$filename"
      else null
    }
    tag "${prefix}"
    label 'perl'

    input:
    tuple val(prefix), val(entrypoint), file(sread1), file(sread2), file(sreads), file(lreads), val(lr_type), file(fast5), file(assembly), val(resfinder_species)
    file(bacannot_db)

    output:
    // Grab all outputs
    path("annotation")
    // Outputs must be linked to each prefix (tag)
    tuple val(prefix), path("annotation/${prefix}.gff")
    tuple val(prefix), path("annotation/${prefix}.gbk")
    tuple val(prefix), path("annotation/${prefix}.fna")
    tuple val(prefix), path("annotation/${prefix}.faa")
    tuple val(prefix), path("annotation/${prefix}.ffn")
    tuple val(prefix), path("annotation/${prefix}.fna"), path("${lreads}"), path("${fast5}")
    tuple val(prefix), path("annotation/${prefix}.fna"), val("${resfinder_species}")
    tuple val(prefix), path("annotation/${prefix}.txt")
    path('prokka_version.txt')

    script:
    kingdom = (params.prokka_kingdom)      ? "--kingdom ${params.prokka_kingdom}"        : ''
    gcode   = (params.prokka_genetic_code) ? "--gcode ${params.prokka_genetic_code}"     : ''
    rnammer = (params.prokka_use_rnammer)  ? "--rnammer"                                 : ''
    pgap    = (params.prokka_skip_pgap)    ? "" : "cp ${bacannot_db}/prokka_db/PGAP_NCBI.hmm \${dbs_dir}/hmm ;"
    """
    # save prokka version
    prokka -v &> prokka_version.txt ;

    # copy additional prokka HMM dbs
    dbs_dir=\$(prokka --listdb 2>&1 >/dev/null |  grep "databases in" | cut -f 4 -d ":" | tr -d " ") ;
    cp ${bacannot_db}/prokka_db/TIGRFAMs_15.0.hmm \${dbs_dir}/hmm ;
    ## get PGAP for prokka if user wants
    ${pgap}                            

    # rebuild prokka databases
    prokka --setupdb ;

    # run prokka
    prokka \\
        $kingdom \\
        $gcode \\
        $rnammer \\
        --outdir annotation \\
        --cpus ${params.threads} \\
        --mincontiglen 200 \\
        --prefix ${prefix} \\
        --genus '' \\
        --species '' \\
        --strain \"${prefix}\" \\
        $assembly
    """
}
