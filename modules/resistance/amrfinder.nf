process amrfinder {
  publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "resistance/AMRFinderPlus/$filename"
  }
  tag "Scanning AMR genes with AMRFinderPlus"
  label 'main'

  input:
  tuple val(prefix), file(proteins)

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

  # Run amrfinder
  CONDA_PREFIX=/opt/conda
  amrfinder -p $proteins --plus -o AMRFinder_complete.tsv --threads ${params.threads} \
  --ident_min \$(echo "scale=2; ${params.blast_resistance_minid}/100" | bc -l ) \
  --coverage_min \$(echo "scale=2; ${params.blast_resistance_mincov}/100" | bc -l ) \
  --name ${prefix} --protein_output ${prefix}_args.faa --database /opt/conda/share/amrfinderplus/data/latest ;
  awk -F '\t' '{ if (\$3 != "") { print } }' AMRFinder_complete.tsv | grep -v "VIRULENCE" > AMRFinder_resistance-only.tsv ;
  """
}
