process AMRFINDER_DB {
    publishDir "${params.output}/amrfinder_db", mode: 'copy', overwrite: "$params.force_update"
    label 'process_ultralow'

    conda "bioconda::ncbi-amrfinderplus=3.11.2"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://ncbi/amr:3.11.2-2022-12-19.1' :
        'ncbi/amr:3.11.2-2022-12-19.1' }"

    output:
    file("*")

    script:
    """
    # download amrfinderplus database
    amrfinder_update -d \$(pwd)
    """
}
