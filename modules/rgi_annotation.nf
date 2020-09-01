process rgi {
  publishDir "${params.outdir}/${prefix}/resistance/RGI", mode: 'copy'
  tag "Resistance genes annotation with RGI"
  label 'main'

  input:
  tuple val(prefix), file(input)

  output:
  // Grab all outputs
  file "*RGI_${prefix}*" optional true
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("Perfect_RGI_${prefix}_hits.txt") optional true
  tuple val(prefix), file("Strict_RGI_${prefix}_hits.txt") optional true
  tuple val(prefix), file("RGI_${prefix}.txt") optional true

  script:
  """
  rgi main --input_sequence $input --output_file ./RGI_${prefix} --input_type protein \
    --alignment_tool DIAMOND --num_threads ${params.threads} --exclude_nudge --clean ;

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
