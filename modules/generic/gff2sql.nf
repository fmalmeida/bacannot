process create_sql {
  publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf(".sqlite") > 0) "sqldb/$filename"
    else "$filename"
  }
  tag "Creating SQL database for the annotation"
  label 'renv'

  input:
    tuple val(prefix), file(gff), file(genes_nt), file(genes_aa), file(genome)

  output:
    file "${prefix}.sqlite"
    file "run_server.sh"

  script:
  """
  # Create SQL db
  gff2sql.R -i $gff -o ${prefix}.sqlite -n $genes_nt -a $genes_aa -f $genome ;

  # Save parser
  cp /work/bscripts/run_server.sh . ;
  """
}
