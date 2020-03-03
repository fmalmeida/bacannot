process mlst {
   publishDir "${params.outdir}/MLST", mode: 'copy'
   container = 'fmalmeida/bacannot:latest'
   tag "Performing MLST analysis"

   input:
     file genome

   output:
     file "${params.prefix}_mlst_analysis.txt" optional true
     file "${params.prefix}_novel_alleles.fasta" optional true

   script:
   """
   source activate MLST ;
   mlst --quiet --novel ${params.prefix}_novel_alleles.fasta $genome > ${params.prefix}_mlst_analysis.txt
   """
}
