process GET_NCBI_PROTEIN {
  label = [ 'misc', 'process_ultralow' ]

  input:
  file(ncbi_accs)

  output:
  path("ncbi_protein.faa")

  script:
  """
  # download and format ncbi protein entries for custom blastp
  for accession in \$(cat ${ncbi_accs}) ; do \\
    esearch -db protein -query "\${accession}" | \\
    efetch -format gp >> ncbi_protein.gbk; \\
  done

  # convert to formatted fasta
  gbk2faa.py ncbi_protein.gbk > ncbi_protein.faa
  """
}
