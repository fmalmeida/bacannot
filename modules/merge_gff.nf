process gff_merge {
  publishDir "${params.outdir}/${prefix}/gffs", mode: 'copy'
  container = 'fmalmeida/bacannot:latest'

  input:
  file gff
  val(prefix)

  output:
  file "${prefix}_merged.gff"

  """
  echo \"##gff-version 3\" > ${prefix}_merged.gff ;
  bedtools sort -i $gff | bedtools merge -d ${params.bedtools_merge_distance} -s -c 2,3,6,7,8,9 -o distinct,distinct,max,distinct,distinct,distinct \
  | awk 'BEGIN { FS = "\t"; OFS="\\t" } sub(",",";",\$9) { print \$1,\$4,\$5,\$2+1,\$3,\$6,\$7,\$8,\$9}' >> ${prefix}_merged.gff
  """
}
