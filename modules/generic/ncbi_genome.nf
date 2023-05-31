process GET_NCBI_GENOME {
    publishDir "${params.output}/closest_genomes", mode: 'copy'
    label = [ 'misc', 'process_low' ]

    input:
    file(ncbi_accs)

    output:
    path("*.fna"), emit: genomes

    when: !params.skip_sourmash

    script:
    """
    # download and format ncbi protein entries for custom blastp
    for biosample in \$(cat ${ncbi_accs}) ; do \\

        acc=\$(esearch -db biosample -query \${biosample} | elink -target assembly | esummary | xtract -pattern DocumentSummary -element Genbank | tr -d '\\n' | tr -d '') && \\
        curl -OJX GET \\
            "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/\${acc}/download?include_annotation_type=GENOME_FASTA&filename=\${acc}.zip" \\
            -H "Accept: application/zip" && \\
        unzip \${acc}.zip && \\
        mv ncbi_dataset/data/*/*.fna . && \\
        rm -rf ncbi_dataset *.zip ;

    done
    """
}
