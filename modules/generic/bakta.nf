process BAKTA {
    publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
      if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
      else if (filename == "annotation") "$filename"
      else null
    }
    tag "${prefix}"
    label = [ 'misc', 'process_high' ]

    input:
    tuple val(prefix), val(entrypoint), file(sread1), file(sread2), file(sreads), file(lreads), val(lr_type), file(fast5), file(assembly), val(resfinder_species)
    file(bakta_db)

    output:
    // Grab all outputs
    file "annotation"
    // Outputs must be linked to each prefix (tag)
    tuple val(prefix), file("annotation/${prefix}.gff3") // annotation in gff format
    tuple val(prefix), file("annotation/${prefix}.gbff") // annotation in gbk format
    tuple val(prefix), file("annotation/${prefix}.fna") // renamed genome
    tuple val(prefix), file("annotation/${prefix}.faa") // gene aa sequences
    tuple val(prefix), file("annotation/${prefix}.ffn") // gene nt sequences
    tuple val(prefix), file("annotation/${prefix}.fna"), file("${lreads}"), file("${fast5}") // For methylation calling
    tuple val(prefix), file("annotation/${prefix}.fna"), val("${resfinder_species}") // For resfinder
    tuple val(prefix), file("annotation/${prefix}.tsv") // bakta stats
    file('bakta_version.txt') // Save bakta version

    script:
    """
    # Save bakta version
    bakta --version &> bakta_version.txt ;

    # Run bakta
    bakta \\
        --output annotation \\
        --threads $task.cpus \\
        --min-contig-length 200 \\
        --prefix ${prefix} \\
        --strain '${prefix}' \\
        --db $bakta_db \\
        $assembly
    """
}
