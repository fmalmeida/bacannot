process rgi {
  publishDir "${params.outdir}/${prefix}/resistance/RGI", mode: 'copy'
  container = 'fmalmeida/bacannot:latest'
  tag "Resistance genes annotation with RGI"

  input:
  file input
  val(prefix)

  output:
  file "Perfect_RGI_${prefix}_hits.txt" optional true
  file "Strict_RGI_${prefix}_hits.txt" optional true
  file "RGI_${prefix}.txt" optional true
  file "*RGI_${prefix}*" optional true

  script:
  """
  source activate RGI ;
  rgi main -i $input -o RGI_${prefix} -t protein -a DIAMOND -n ${params.threads} --exclude_nudge --clean ;

  ## Parse perfect hits
  sed 's/ # /#/g' RGI_${prefix}.txt \
  | awk 'BEGIN { FS = "\t"; OFS="\\t" } ; { print \$2,\$3,\$4,\$5,\$6,\$9,\$11,\$15,\$16,\$17 }' \
  | awk '{ if (\$5 == "Perfect") print }' > Perfect_RGI_${prefix}_hits.txt

  ## Parse strict hits
  sed 's/ # /#/g' RGI_${prefix}.txt \
  | awk 'BEGIN { FS = "\t"; OFS="\\t" } ; { print \$2,\$3,\$4,\$5,\$6,\$9,\$11,\$15,\$16,\$17 }' \
  | awk '{ if (\$5 == "Strict") print }' > Strict_RGI_${prefix}_hits.txt
  """
}
