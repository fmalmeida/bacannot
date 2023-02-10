process KOFAMSCAN {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "$filename"
  }
  tag "${prefix}"
  label = [ 'process_high' ]

  conda "bioconda::kofamscan=1.3.0"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/kofamscan:1.3.0--hdfd78af_2' :
        'quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2' }"

  input:
  tuple val(prefix), file('proteins.faa')
  file(bacannot_db)

  output:
  // Grab all outputs
  tuple val(prefix), path("KOfamscan"), emit: all
  tuple val(prefix), path("KOfamscan/${prefix}_ko_forKEGGMapper.txt"), emit: results

  script:
  """
  # Get kofamscan version
  exec_annotation -v | sed "s/exec_annotation/kofamscan/" > kofamscan_version.txt

  # Create dir for results
  mkdir KOfamscan ;

  # Run kofamscan with detailed output
  exec_annotation \\
      -p ${bacannot_db}/kofamscan_db/profiles/prokaryote.hal \\
      -k ${bacannot_db}/kofamscan_db/ko_list \\
      -o KOfamscan/${prefix}_ko_detailed.txt \\
      --keep-tabular \\
      --cpu=$task.cpus \\
      proteins.faa ;

  # Re-run kofamscan with mapper-output
  exec_annotation \\
      -p ${bacannot_db}/kofamscan_db/profiles/prokaryote.hal \\
      -k ${bacannot_db}/kofamscan_db/ko_list \\
      -o KOfamscan/${prefix}_ko_forKEGGMapper.txt \\
      --reannotate \\
      --cpu=$task.cpus \\
      -f mapper \\
      proteins.faa ;
  """
}
