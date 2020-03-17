process prokka {
    publishDir "${params.outdir}/${prefix}", mode: 'copy'
    container = 'fmalmeida/bacannot:latest'
    tag "Executing generic gene annotation with Prokka"

    input:
    file input
    val(prefix)

    output:
    file "prokka/${prefix}.*" // needed to take all output into the output dir
    file "prokka/${prefix}.gff" // annotation in gff format
    file "prokka/${prefix}.gbk" // annotation in gbk format
    file "prokka/${prefix}.fna" // renamed genome
    file "prokka/${prefix}.faa" // gene aa sequences
    file "prokka/${prefix}.ffn" // gene nt sequences
    val(prefix) // Save genome filename

    script:
    kingdom = (params.prokka_kingdom) ? "--kingdom ${params.prokka_kingdom}" : ''
    gcode = (params.prokka_genetic_code) ? "--gcode ${params.prokka_genetic_code}" : ''
    rnammer = (params.prokka_use_rnammer) ? "--rnammer" : ''
    genus = (params.prokka_genus) ? "--genus ${params.prokka_genus} --usegenus" : ''
    """
    source activate PROKKA ;
    prokka $kingdom $gcode $rnammer --outdir prokka --cpus ${params.threads} --centre ${params.prokka_center} \
    --mincontiglen 200 $genus --prefix ${prefix} $input
    """
}
