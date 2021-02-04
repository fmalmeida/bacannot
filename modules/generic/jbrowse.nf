process jbrowse {
  publishDir "${params.outdir}/${prefix}/jbrowse", mode: 'copy'
  label 'jbrowse'
  tag "Creating the genome browser with JBrowse"

  input:
  tuple val(prefix), file(gff), file(draft), file("prokka_gff"), file(mlst), file(barrnap),
        file(gc_bedGraph), file(gc_chrSizes), file(kofamscan), file(vfdb),
        file(victors), file(amrfinder), file(rgi), file(iceberg), file(phast),
        file(phigaro), file(genomic_islands), file("methylation"), file("chr.sizes"),
        file(phispy_tsv), file(resfinder_gff)

  output:
  file "*"

  """
  # Get JBrowse Files in working directory
  cp -R /work/jbrowse/* . ;
  cp /work/bscripts/run_jbrowse.sh . ;

  # Render genome browser
  ./run_jbrowse.sh -p $prefix -g $draft -b $gc_bedGraph -s $gc_chrSizes -f $gff -r $barrnap -B $phigaro \
  -P $phispy_tsv -G $genomic_islands -m methylation -S chr.sizes -R $resfinder_gff
  """
}
