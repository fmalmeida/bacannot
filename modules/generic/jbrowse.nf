process JBROWSE {
  publishDir "${params.output}/${prefix}/jbrowse", mode: 'copy'
  label = [ 'jbrowse', 'process_low' ]
  tag "${prefix}"

  input:
  tuple val(prefix), file(merged_gff), file(draft), file("prokka_gff"), file(barrnap), file(gc_bedGraph), file(gc_chrSizes), file(resfinder_gff), file(phigaro), file(genomic_islands), file("methylation"), file("chr.sizes"), file(phispy_tsv), file(digIS_gff), file(antiSMASH), file(custom_annotations), file(integron_finder)

  output:
  path "*", emit: results

  script:
  """
  # Get JBrowse Files in working directory
  cp -R /work/jbrowse/* . ;

  # Render genome browser
  run_jbrowse.sh \\
    -p $prefix \\
    -g $draft \\
    -b $gc_bedGraph \\
    -s $gc_chrSizes \\
    -f $merged_gff \\
    -r $barrnap \\
    -B $phigaro \\
    -P $phispy_tsv \\
    -G $genomic_islands \\
    -m methylation \\
    -S chr.sizes \\
    -R $resfinder_gff \\
    -d $digIS_gff \\
    -A $antiSMASH \\
    -i $integron_finder
  """
}
