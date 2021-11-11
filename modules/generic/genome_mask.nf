process masking_genome {
  publishDir "${params.outdir}/${prefix}", mode: 'copy',
  saveAs: {filename ->
  //This line saves the files with specific sufixes in specific folders
  if (filename.indexOf(".gff") > 0 ) "gffs/$filename"
  else if (filename.indexOf(".fasta") > 0 ) "masked_genome/$filename"
  else if (filename.indexOf(".txt") > 0 ) "gffs/$filename"
}
  tag "${prefix}"
  label 'main'

  input:
  file input
  file 'gff'
  val(prefix)

  output:
  file "${prefix}_clear.gff" // Annotation in GFF format without the genome sequences
  file "${prefix}_masked_genome.fasta" // Masked genome in FASTA file
  file "readme.txt"

  """
  grep "ID=" gff | awk '{ print \$1 "\t" \$4 "\t" \$5 }' > cds_prokka.bed ;
  grep "ID=" gff > ${prefix}_clear.gff ;
  maskFastaFromBed -fi $input -fo ${prefix}_masked_genome.fasta -bed cds_prokka.bed ;
  echo -e \"Understanding the *_clear.gff file.\\n${prefix}_clear.gff is the same GFF file from Prokka output. However, this file does \
  not contain any genomic sequences as Prokka output does.\" > readme.txt
  """
}
