process PHISPY {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "PhiSpy") "prophages/$filename"
    else null
  }
  tag "${prefix}"
  label = [ 'python', 'process_medium' ]

  input:
  tuple val(prefix), file(input)

  output:
  tuple val(prefix), path("PhiSpy")             , emit: all
  tuple val(prefix), path("PhiSpy/prophage.tsv"), emit: results
  tuple val(prefix), path("phispy_version.txt") , emit: version

  script:
  """
  # get tool version
  PhiSpy.py -v > phispy_version.txt ;

  # run phispy
  PhiSpy.py \\
      -o PhiSpy \\
      $input \\
      --color \\
      --output_choice 127 \\
      --threads $task.cpus
  """
}
