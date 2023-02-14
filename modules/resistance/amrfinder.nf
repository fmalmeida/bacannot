process AMRFINDER {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "resistance/AMRFinderPlus/$filename"
  }
  tag "${prefix}"
  label = [ 'process_medium' ]

  conda "bioconda::ncbi-amrfinderplus=3.11.2"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://ncbi/amr:3.11.2-2022-12-19.1' :
        'ncbi/amr:3.11.2-2022-12-19.1' }"

  input:
  tuple val(prefix), file(proteins)
  file(bacannot_db)

  output:
  tuple val(prefix), path("AMRFinder_resistance-only.tsv"), emit: resistance_results
  tuple val(prefix), path("AMRFinder_complete.tsv")       , emit: complete_results
  tuple val(prefix), path("*")                            , emit: all
  path("${prefix}_args.faa")                              , emit: proteins
  path("amrfinder_version.txt")                           , emit: version

  script:
  resistance_minid  = params.blast_resistance_minid / 100.00
  resistance_mincov = params.blast_resistance_mincov / 100.00
  """
  # Get tool version
  amrfinder --version > amrfinder_version.txt ;

  # run amrfinder
  amrfinder \\
      -p $proteins \\
      --plus \\
      -o AMRFinder_complete.tsv \\
      --threads $task.cpus \\
      --ident_min ${resistance_minid} \\
      --coverage_min ${resistance_mincov} \\
      --name ${prefix} \\
      --protein_output ${prefix}_args.faa \\
      --database ${bacannot_db}/amrfinder_db/latest
  
  # filter results
  awk \\
      -F '\t' \\
      '{ if (\$3 != "") { print } }' \\
      AMRFinder_complete.tsv | \\
      grep -v "VIRULENCE" > AMRFinder_resistance-only.tsv ;
  """
}
