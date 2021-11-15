process jbrowse {
  publishDir "${params.output}/${prefix}/jbrowse", mode: 'copy'
  label 'jbrowse'
  tag "${prefix}"

  input:
  tuple val(prefix), file(gff), file(draft), file("prokka_gff"), file(mlst), file(barrnap),
        file(gc_bedGraph), file(gc_chrSizes), file(kofamscan), file(vfdb),
        file(victors), file(amrfinder), file(resfinder_gff), file(rgi), file(iceberg), file(phast),
        file(phigaro), file(genomic_islands), file("methylation"), file("chr.sizes"),
        file(phispy_tsv), file("digIS.gff"), file("antiSMASH")

  output:
  file "*"

  script:
  """
  # Get JBrowse Files in working directory
  cp -R /work/jbrowse/* . ;
  cp /work/bscripts/run_jbrowse.sh . ;
  chmod a+x run_jbrowse.sh ;

  # Render genome browser
  ./run_jbrowse.sh -p $prefix -g $draft -b $gc_bedGraph -s $gc_chrSizes -f $gff -r $barrnap -B $phigaro \
  -P $phispy_tsv -G $genomic_islands -m methylation -S chr.sizes -R $resfinder_gff -d digIS.gff
  """
}
