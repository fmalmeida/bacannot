process GET_NCBI_GENOME {
    publishDir "${params.output}/sourmash_all/closest_genomes", mode: 'copy'
    label = [ 'misc', 'process_low' ]
    maxForks 3

    input:
    val biosample_acc

    output:
    path("*.fna"), emit: genomes

    when: !params.skip_sourmash

    script:
    """
    # download genomes from ncbi
    export acc=\$(esearch -db biosample -query '${biosample_acc}' | elink -target assembly | esummary | xtract -pattern DocumentSummary -element Genbank | tr -d '\\n' | tr -d '')
    echo ${biosample_acc} -- \${acc}
    curl -OJX GET \\
        --retry 5 \\
        --retry-delay 15 \\
        "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/\${acc}/download?include_annotation_type=GENOME_FASTA&filename=\${acc}.zip" \\
        -H "Accept: application/zip"

    # unzip files
    for file in *.zip ; do
        rm -rf ncbi_dataset *.md *.txt && \\
        unzip -u \$file && \\
        mv ncbi_dataset/data/*/*.fna . && \\
        rm -rf ncbi_dataset *.md *.txt
    done

    # rename
    for file in *.fna ; do
        name=\$( echo \$file | cut -d '_' -f 1,2 ) ;
        mv \$file \${name}.fna
    done
    """
}
