process masking_genome {
  container 'fmalmeida/bacannot:latest'
  tag "Masking genome with bedtools"

  input:
  file input
  file 'gff'

  output:
  file "${params.prefix}_clear.gff" // Annotation in GFF format without the genome sequences
  file "${params.prefix}_masked_genome.fasta" // Masked genome in FASTA file

  """
  grep "ID=" gff | awk '{ print \$1 "\t" \$4 "\t" \$5 }' > cds_prokka.bed ;
  grep "ID=" gff > ${params.prefix}_clear.gff ;
  maskFastaFromBed -fi $input -fo ${params.prefix}_masked_genome.fasta -bed cds_prokka.bed
  """
}
