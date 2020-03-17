process barrnap {
   publishDir "${params.outdir}/${prefix}/rRNA", mode: 'copy'
   container = 'fmalmeida/bacannot:latest'
   tag "Predicting rRNA sequences with barrnap pipeline from T. Seeman"

   input:
   file input
   val(prefix)

   output:
   file "${prefix}_rRNA.gff"
   file "${prefix}_rRNA.fa" optional true

   script:
   """
   barrnap -o ${prefix}_rRNA.fa < $input > ${prefix}_rRNA.gff
   """
}
