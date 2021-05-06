process mlst {
   publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
     if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
     else "MLST/$filename"
   }
   tag "Performing MLST analysis with mlst pipeline from T. Seeman"
   label 'main'

   input:
   tuple val(prefix), file(genome)

   output:
   tuple val(prefix), file("${prefix}_mlst_analysis.txt") optional true
   tuple val(prefix), file("${prefix}_novel_alleles.fasta") optional true
   file('mlst_version.txt')

   script:
   """
   # activate env
   source activate PERL_env ;

   # Save mlst tool version
   mlst --version > mlst_version.txt ;

   # Run mlst
   mlst --quiet --novel ${prefix}_novel_alleles.fasta $genome > ${prefix}_mlst_analysis.txt
   """
}
