process prokka {
    publishDir "${params.outdir}", mode: 'copy'
    container = 'fmalmeida/bacannot:latest'

    input:
    file input
    val threads

    output:
    file "prokka/${params.prefix}.*"
    file "prokka/${params.prefix}.gff" // annotation in gff format
    file "prokka/${params.prefix}.gbk" // annotation in gbk format
    file "prokka/${params.prefix}.fna" // renamed genome
    file "prokka/${params.prefix}.faa" // gene aa sequences
    file "prokka/${params.prefix}.ffn" // gene nt sequences

    script:
    kingdom = (params.prokka_kingdom) ? "--kingdom ${params.prokka_kingdom}" : ''
    gcode = (params.prokka_genetic_code) ? "--gcode ${params.prokka_genetic_code}" : ''
    rnammer = (params.prokka_use_rnammer) ? "--rnammer" : ''
    genus = (params.prokka_genus) ? "--genus ${params.prokka_genus} --usegenus" : ''
    """
    source activate PROKKA ;
    prokka $kingdom $gcode $rnammer --outdir prokka --cpus $threads --centre ${params.prokka_center} \
    --mincontiglen 200 $genus --prefix ${params.prefix} $input
    """
}
