process PHISPY {
  publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "PhiSpy") "prophages/$filename"
    else null
  }
  tag "${prefix}"

  input:
  tuple val(prefix), file(input)

  output:
  tuple val(prefix), path("PhiSpy"), emit: all
  tuple val(prefix), path("PhiSpy/prophage.tsv"), emit: tsv
  tuple val(prefix), path("phispy_version.txt") , emit: txt

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
      --threads ${params.threads}
  """
}
