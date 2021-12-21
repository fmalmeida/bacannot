process KOFAMSCAN {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "$filename"
  }
  tag "${prefix}"

  input:
  tuple val(prefix), file('proteins.faa')
  file(bacannot_db)

  output:
  // Grab all outputs
  file("KOfamscan")
  tuple val(prefix), file("KOfamscan/${prefix}_ko_forKEGGMapper.txt")

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
      --cpu=${params.threads} \\
      proteins.faa ;

  # Re-run kofamscan with mapper-output
  exec_annotation \\
      -p ${bacannot_db}/kofamscan_db/profiles/prokaryote.hal \\
      -k ${bacannot_db}/kofamscan_db/ko_list \\
      -o KOfamscan/${prefix}_ko_forKEGGMapper.txt \\
      --reannotate \\
      --cpu=${params.threads} \\
      -f mapper \\
      proteins.faa ;
  """
}
