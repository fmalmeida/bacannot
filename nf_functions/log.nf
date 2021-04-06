def logMessage() {
  log.info "================================================================="
  log.info " Container-based, fmalmeida/bacannot, Genome Annotation Pipeline "
  log.info "================================================================="
  def summary = [:]
  if (params.genome) { summary['Input genomes'] = params.genome }
  summary['Output dir']   = "${params.outdir}"
  summary['Threads'] = params.threads
  if (params.skip_virulence_search == false) {
  summary['Blast % ID - Virulence Genes'] = params.blast_virulence_minid
  summary['Blast query coverage - Virulence Genes'] = params.blast_virulence_mincov
  }
  if (params.skip_resistance_search == false) {
  summary['Blast % ID - AMR Genes'] = params.blast_resistance_minid
  summary['Blast query coverage - AMR Genes'] = params.blast_resistance_mincov
  }
  if (params.skip_iceberg_search == false | params.skip_prophage_search == false) {
  summary['Blast % ID - ICEs or Phages'] = params.blast_MGEs_minid
  summary['Blast query coverage - ICEs or Phages'] = params.blast_MGEs_mincov
  }
  if (params.skip_plasmid_search == false) {
  summary['Blast % ID - Plasmids'] = params.plasmids_minid
  summary['Blast query coverage - Plasmids'] = params.plasmids_mincov
  }
  if(workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Current home']   = "$HOME"
  summary['Current user']   = "$USER"
  summary['Current path']   = "$PWD"
  summary['Configuration file'] = workflow.configFiles[0]
  log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
  log.info "=============================================================="
}
