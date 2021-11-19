process PROKKA {
    publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
      if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
      else if (filename == "annotation") "$filename"
      else null
    }
    tag "${prefix}"

    input:
    tuple val(prefix), val(entrypoint), file(sread1), file(sread2), file(sreads), file(lreads), val(lr_type), file(fast5), file(assembly), val(resfinder_species)
    file(bacannot_db)

    output:
    // Grab all outputs
    file "annotation"
    // Outputs must be linked to each prefix (tag)
    tuple val(prefix), path("annotation/${prefix}.gff"), emit: gff
    tuple val(prefix), path("annotation/${prefix}.gbk"), emit: gbk
    tuple val(prefix), path("annotation/${prefix}.fna"), emit: renamedGenome
    tuple val(prefix), path("annotation/${prefix}.faa"), emit: genesAA
    tuple val(prefix), path("annotation/${prefix}.ffn"), emit: genesNT
    tuple val(prefix), path("annotation/${prefix}.fna"), path("${lreads}"), path("${fast5}"), emit: fast5
    tuple val(prefix), path("annotation/${prefix}.fna"), val("${resfinder_species}"), emit: resfinder
    tuple val(prefix), path("annotation/${prefix}.txt"), emit: stats
    path('prokka_version.txt'), emit: version

    script:
    kingdom = (params.prokka_kingdom)      ? "--kingdom ${params.prokka_kingdom}"        : ''
    gcode   = (params.prokka_genetic_code) ? "--gcode ${params.prokka_genetic_code}"     : ''
    rnammer = (params.prokka_use_rnammer)  ? "--rnammer"                                 : ''
    """
    # save prokka version
    prokka -v &> prokka_version.txt ;

    # rebuild prokka dbs with downloaded HMM
    cp ${bacannot_db}/prokka_db/* \$(find / -name "hmm" -type d) ;
    prokka --setupdb ;

    # Run prokka
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
