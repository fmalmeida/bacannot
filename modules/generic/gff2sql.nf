process create_sql {
  publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf(".sqlite") > 0) "sqldb/$filename"
    else "$filename"
  }
  tag "Creating SQL database for the annotation"
  label 'renv'

  input:
    tuple val(prefix), file(gff), file(genes_nt), file(genes_aa), file(genome), file("digIS.gff"), file("digIS.fa"), file("digIS.faa")

  output:
    file "${prefix}.sqlite"
    file "run_server.sh"

  script:
  """
  if [ -s digIS.fa ] ;
  then

    # concatenate files
    cat $gff digIS.gff | bedtools sort > input.gff
    cat $genes_nt digIS.fa  > input.fa
    cat $genes_aa digIS.faa > input.faa

    # Create SQL db
    gff2sql.R -i input.gff -o ${prefix}.sqlite -n input.fa -a input.faa -f $genome &> gff2sql.log;

  else

    # Create SQL db
    gff2sql.R -i $gff -o ${prefix}.sqlite -n $genes_nt -a $genes_aa -f $genome &> gff2sql.log;

  fi

  # Save results with better name
  mv /work/${prefix}.sqlite . ;

  # Save parser
  cp /work/bscripts/run_server.sh . ;
  """
}
