process PHIGARO {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "prophages/phigaro/$filename"
  }
  tag "${prefix}"
  label = [ 'process_medium' ]

  conda "bioconda::phigaro=2.3.0"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/phigaro:2.3.0--pyh7b7c402_0' :
        'quay.io/biocontainers/phigaro:2.3.0--pyh7b7c402_0' }"

  input:
  tuple val(prefix), file("assembly.fasta")
  file(bacannot_db)

  output:
  tuple val(prefix), path("${prefix}_phigaro.tsv") , emit: tsv
  tuple val(prefix), path("${prefix}_phigaro.bed") , emit: bed
  tuple val(prefix), path("${prefix}_phigaro.html"), emit: html optional true
  tuple val(prefix), path("*")                     , emit: all
  path('phigaro_version.txt')                      , emit: version

  script:
  """  
  # get tool version
  phigaro -V > phigaro_version.txt ;

  # create new config to properly load database
  cp \$(which config.yml) ./custom_config.yml ;
  sed -i "s|CHANGE_PVOG|${bacannot_db}/phigaro_db/allpvoghmms|" ./custom_config.yml ;
  HMM_BIN=\$(which hmmsearch) ;
  sed -i "s|CHANGE_HMMSEARCH|\$HMM_BIN|" ./custom_config.yml ;
  PRODIGAL_BIN=\$(which prodigal) ;
  sed -i "s|CHANGE_PRODIGAL|\$PRODIGAL_BIN|" ./custom_config.yml ;

  # run phigaro
  phigaro \\
      -f assembly.fasta \\
      --config ./custom_config.yml \\
      -t $task.cpus \\
      -e html tsv \\
      -o out.phg \\
      --delete-shorts \\
      -p \\
      --not-open ;

  # change names
  [ ! -s out.phg/assembly.phigaro.tsv  ] || cp out.phg/assembly.phigaro.tsv ${prefix}_phigaro.tsv ;
  [ ! -s out.phg/assembly.phigaro.html ] || cp out.phg/assembly.phigaro.html ${prefix}_phigaro.html ;

  # create BED
  touch ${prefix}_phigaro.bed ;
  [ ! -s out.phg/assembly.phigaro.html ] || grep -v "taxonomy" ${prefix}_phigaro.tsv | \\
      awk 'BEGIN { FS = "\\t"; OFS="\\t" } { print \$1,\$2,\$3 }' > ${prefix}_phigaro.bed
  """
}
