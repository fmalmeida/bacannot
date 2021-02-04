process call_methylation {
  publishDir "${params.outdir}/${prefix}", mode: 'copy', saveAs: { filename ->
    if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
    else "methylations/$filename"
  }
  tag "5mC (CpG) Methylated sites are being calculated with Nanopolish"
  label 'main'

  input:
  tuple val(prefix), file(draft)
  file(fast5)
  file(reads)

  output:
  // Grab all outputs
  file "*_calls.tsv" optional true
  file "*_frequency.tsv" optional true
  tuple val(prefix), file("methylation_frequency.bedGraph") optional true
  tuple val(prefix), file("chr.sizes") optional true
  file('nanopolish_version.txt')

  script:
  fast5_dir = fast5.getName()
  """
  # Get tool version
  nanopolish --version > nanopolish_version.txt ;

  # Index Our Fast5 Data
  nanopolish index -d ${fast5_dir} ${reads} ;

  # Map Our Indexed Reads to Our Genome
  minimap2 -a -x map-ont ${draft} ${reads} | samtools sort -T tmp -o reads_output.sorted.bam ;
  samtools index reads_output.sorted.bam ;

  # Call Methylation
  nanopolish call-methylation -r ${reads} -b reads_output.sorted.bam -g ${draft} -t ${params.threads} > methylation_call.tsv ;

  # Calculate Methylation Frequencies
  /work/nanopolish/scripts/calculate_methylation_frequency.py methylation_call.tsv > methylation_frequency.tsv ;

  # Transform These TSV files into bedGraph
  [ ! -s methylation_frequency.tsv ] || grep -v "start" methylation_frequency.tsv | \
  awk '{ print \$1 "\t" \$2 "\t" \$3 "\t" \$7 }' > methylation_frequency.bedGraph ;

  # Create Contig Sizes File
  seqtk comp ${draft} | awk '{ print \$1 "\t" \$2 }' > chr.sizes
  """
}
