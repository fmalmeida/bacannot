process BARRNAP {
   publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
     if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
     else "rRNA/$filename"
   }
   tag "${prefix}"
   label 'perl'

   input:
   tuple val(prefix), file(genome)

   output:
   tuple val(prefix), path("${prefix}_rRNA.gff")
   tuple val(prefix), path("${prefix}_rRNA.fa")
   path('barrnap_version.txt')

   script:
   """
   # save barrnap tool version
   barrnap --version &> barrnap_version.txt ;

   # run barrnap
   barrnap -o ${prefix}_rRNA.fa < $genome > ${prefix}_rRNA.gff
   """
}
