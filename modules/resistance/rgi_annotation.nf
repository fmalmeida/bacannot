process CARD_RGI {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "Parsed_RGI_${prefix}_hits.txt") null
    else "resistance/RGI/$filename"
  }
  tag "${prefix}"
  label = [ 'python', 'process_medium' ]

  input:
  tuple val(prefix), path(input)
  path(bacannot_db)

  output:
  path "*RGI_${prefix}*"                                  , emit: all         optional true 
  tuple val(prefix), path("Parsed_RGI_${prefix}_hits.txt"), emit: parsed_hits optional true
  tuple val(prefix), path("RGI_${prefix}.txt")            , emit: raw_hits    optional true 
  tuple val(prefix), path("heatmap/RGI*heatmap*.png")     , emit: heatmap_png optional true
  path("heatmap")                                         , emit: heatmap     optional true 
  path("*_version.txt")                                   , emit: version

  script:
  """
  # activate env
  source activate rgi
  
  # load database
  rgi load --card_json ${bacannot_db}/card_db/card.json --local

  # get tool version
  rgi main --version > rgi_version.txt ;
  rgi database --version --local > card_db_version.txt ;

  # execute RGI
  rgi main \\
      --local \\
      --input_sequence $input \\
      --output_file RGI_${prefix}_unfiltered \\
      --input_type protein \\
      --num_threads $task.cpus \\
      --clean ;

  ## filtering by identity
  awk '
      BEGIN { FS = "\\t"; OFS="\\t" } 
      { if (\$10 >= ${params.blast_resistance_minid}) print }
      ' RGI_${prefix}_unfiltered.txt > ./RGI_${prefix}.txt

  ## parse RGI results for reports
  cat RGI_${prefix}.txt | \\
      tail -n +2 | \\
      awk '
          BEGIN { FS = "\\t"; OFS="\\t" }
          { split(\$1,a," "); print a[1],\$6,\$9,\$11,\$15,\$16,\$17 }
          ' > Parsed_RGI_${prefix}_hits.txt

  # draw heatmap for single sample
  if [ \$(wc -l RGI_${prefix}.txt | cut -d " " -f 1) -gt 1 ]
  then
    mkdir -p heatmap ;
    cp RGI_${prefix}_unfiltered.json heatmap/${prefix}.json ;
    rgi heatmap --input ./heatmap -cat drug_class -d text ;
    rm heatmap/${prefix}.json ;
    mv RGI*heatmap* heatmap ;
  fi
  """
}
