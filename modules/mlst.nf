process mlst {
   publishDir "${params.outdir}/${prefix}/MLST", mode: 'copy'
   container = 'fmalmeida/bacannot:latest'
   tag "Performing MLST analysis with mlst pipeline from T. Seeman"

   input:
     tuple val(prefix), file(genome)

   output:
     tuple val(prefix), file("${prefix}_mlst_analysis.txt") optional true
     tuple val(prefix), file("${prefix}_novel_alleles.fasta") optional true

   script:
   """
   source activate MLST ;
   mlst --quiet --novel ${prefix}_novel_alleles.fasta $genome > ${prefix}_mlst_analysis.txt
   """
}
