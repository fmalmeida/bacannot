process mlst {
   publishDir "${params.outdir}/${prefix}/MLST", mode: 'copy'
   container = 'fmalmeida/bacannot:dev'
   tag "Performing MLST analysis with mlst pipeline from T. Seeman"

   input:
     file(genome)

   output:
     tuple val(prefix), file("${prefix}_mlst_analysis.txt") optional true
     tuple val(prefix), file("${prefix}_novel_alleles.fasta") optional true

   script:
   prefix = "${genome.baseName}"
   """
   source activate MLST ;
   mlst --quiet --novel ${prefix}_novel_alleles.fasta $genome > ${prefix}_mlst_analysis.txt
   """
}
