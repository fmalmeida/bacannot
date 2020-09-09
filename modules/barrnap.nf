process barrnap {
   publishDir "${params.outdir}/${prefix}/rRNA", mode: 'copy'
   tag "Predicting rRNA sequences with barrnap pipeline from T. Seeman"
   label 'main'

   input:
   tuple val(prefix), file(genome)

   output:
   tuple val(prefix), file("${prefix}_rRNA.gff")
   tuple val(prefix), file("${prefix}_rRNA.fa") optional true

   script:
   """
   barrnap -o ${prefix}_rRNA.fa < $genome > ${prefix}_rRNA.gff
   """
}
