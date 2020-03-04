process barrnap {
   publishDir "${params.outdir}/rRNA", mode: 'copy'
   container = 'fmalmeida/bacannot:latest'

   input:
   file input

   output:
   file "${params.prefix}_rRNA.gff"
   file "${params.prefix}_rRNA.fa" optional true

   script:
   """
   barrnap -o ${params.prefix}_rRNA.fa < $input > ${params.prefix}_rRNA.gff
   """
}
