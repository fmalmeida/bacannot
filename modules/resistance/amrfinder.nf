process AMRFINDER {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "resistance/AMRFinderPlus/$filename"
  }
  tag "${prefix}"
  label = [ 'misc', 'process_medium' ]

  input:
  tuple val(prefix), file(proteins)
  file(bacannot_db)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("AMRFinder_resistance-only.tsv")
  tuple val(prefix), file("AMRFinder_complete.tsv")
  file("${prefix}_args.faa")
  file("amrfinder_version.txt")

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
      --threads ${params.threads} \\
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
