process barrnap {
   publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
     if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
     else "rRNA/$filename"
   }
   tag "Predicting rRNA sequences with barrnap pipeline from T. Seeman"
   label 'main'

   input:
   tuple val(prefix), file(genome)

   output:
   tuple val(prefix), file("${prefix}_rRNA.gff")
   tuple val(prefix), file("${prefix}_rRNA.fa") optional true
   file('barrnap_version.txt')

   script:
   """
   # Save barrnap tool version
   barrnap --version> barrnap_version.txt ;

   # Run barrnap
   barrnap -o ${prefix}_rRNA.fa < $genome > ${prefix}_rRNA.gff
   """
}
