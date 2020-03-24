process prokka {
    publishDir "${params.outdir}/${prefix}", mode: 'copy'
    container = 'fmalmeida/bacannot:latest'
    tag "Executing generic gene annotation with Prokka"

    input:
    tuple val(prefix), file(input)

    output:
    // Grab all outputs
    file "prokka/${prefix}.*"
    // Outputs must be linked to each prefix (tag)
    tuple val(prefix), file("prokka/${prefix}.gff") // annotation in gff format
    tuple val(prefix), file("prokka/${prefix}.gbk") // annotation in gbk format
    tuple val(prefix), file("prokka/${prefix}.fna") // renamed genome
    tuple val(prefix), file("prokka/${prefix}.faa") // gene aa sequences
    tuple val(prefix), file("prokka/${prefix}.ffn") // gene nt sequences

    script:
    kingdom = (params.prokka_kingdom)      ? "--kingdom ${params.prokka_kingdom}"        : ''
    gcode   = (params.prokka_genetic_code) ? "--gcode ${params.prokka_genetic_code}"     : ''
    rnammer = (params.prokka_use_rnammer)  ? "--rnammer"                                 : ''
    genus   = (params.prokka_genus)        ? "--genus ${params.prokka_genus} --usegenus" : ''
    """
    source activate PROKKA ;
    prokka $kingdom $gcode $rnammer --outdir prokka --cpus ${params.threads} --centre ${params.prokka_center} \
    --mincontiglen 200 $genus --prefix ${prefix} $input
    """
}
