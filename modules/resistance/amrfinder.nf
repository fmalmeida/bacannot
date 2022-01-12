process AMRFINDER {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "resistance/AMRFinderPlus/$filename"
  }
  tag "${prefix}"

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
  """
  # Get tool version
  amrfinder --version > amrfinder_version.txt ;

  # run amrfinder
  amrfinder \\
      -p $proteins \\
      --plus \\
      -o AMRFinder_complete.tsv \\
      --threads ${params.threads} \\
      --ident_min \$(echo print ${params.blast_resistance_minid}/100" | perl ) \\
      --coverage_min \$(echo print ${params.blast_resistance_mincov}/100" | perl ) \\
      --name ${prefix} \\
      --protein_output ${prefix}_args.faa \\
      --database ${bacannot_db}/amrfinder_db
  
  # filter results
  awk -F '\t' '{ if (\$3 != "") { print } }' AMRFinder_complete.tsv | \\
      grep -v "VIRULENCE" > AMRFinder_resistance-only.tsv ;
  """
}
