process phispy {
  publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else if (filename == "PhiSpy") "prophages/$filename"
    else null
  }
  tag "Detecting prophages with PhiSpy"
  label 'main'

  input:
  tuple val(prefix), file(input)

  output:
  tuple val(prefix), file("PhiSpy")
  tuple val(prefix), file("PhiSpy/prophage.tsv")
  tuple val(prefix), file("phispy_version.txt")

  script:
  """
  # Activate conda environment
  source activate phispy ;

  # Get tool version
  PhiSpy.py -v > phispy_version.txt ;

  # Run phispy
  PhiSpy.py -o PhiSpy $input --color --output_choice 127 --threads ${params.threads}
  """
}
