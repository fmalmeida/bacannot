process card_rgi {
  publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "resistance/RGI/$filename"
  }
  tag "Scanning AMR genes with RGI"
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
  tuple val(prefix), file("RGI*heatmap*.png") optional true
  file("heatmap")
  file("*_version.txt")

  script:
  """
  # Activate environment
  source activate RGI ;

  # Get tool version
  rgi main --version > rgi_version.txt ;
  rgi database --version > card_db_version.txt ;

  # Execute RGI
  rgi main --input_sequence $input --output_file RGI_${prefix}_unfiltered --input_type protein \
  --num_threads ${params.threads} --exclude_nudge --clean ;

  ## Filtering by identity
  awk 'BEGIN { FS = "\\t"; OFS="\\t" } { if (\$10 >= ${params.blast_resistance_minid}) print }' RGI_${prefix}_unfiltered.txt > ./RGI_${prefix}.txt

  ## Parse perfect hits
  cat RGI_${prefix}.txt  | \
  awk 'BEGIN { FS = "\\t"; OFS="\\t" } ; { split(\$1,a," "); print a[1],\$6,\$9,\$11,\$15,\$16,\$17 }' | \
  awk '{ if (\$2 == "Perfect") print }'  > Perfect_RGI_${prefix}_hits.txt

  ## Parse strict hits
  cat RGI_${prefix}.txt  | \
  awk 'BEGIN { FS = "\\t"; OFS="\\t" } ; { split(\$1,a," "); print a[1],\$6,\$9,\$11,\$15,\$16,\$17 }' | \
  awk '{ if (\$2 == "Strict") print }'  > Strict_RGI_${prefix}_hits.txt

  # Draw heatmap for single sample
  mkdir -p heatmap ;
  cp RGI_${prefix}_unfiltered.json heatmap/${prefix}.json ;
  rgi heatmap --input ./heatmap -cat drug_class -d text ;
  rm heatmap/${prefix}.json
  mv RGI*heatmap* heatmap ;
  """
}
