process prokka {
    publishDir "${params.outdir}/${prefix}", mode: 'copy'
    container = 'fmalmeida/bacannot:dev'
    tag "Executing generic annotation with Prokka"

    input:
    file(input)

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
    prefix  = "${input.baseName}"
    kingdom = (params.prokka_kingdom)      ? "--kingdom ${params.prokka_kingdom}"        : ''
    gcode   = (params.prokka_genetic_code) ? "--gcode ${params.prokka_genetic_code}"     : ''
    rnammer = (params.prokka_use_rnammer)  ? "--rnammer"                                 : ''
    """
    source activate PROKKA ;
    prokka $kingdom $gcode $rnammer --outdir prokka \
    --cpus ${params.threads} --mincontiglen 200 --prefix ${prefix} $input
    """
}
