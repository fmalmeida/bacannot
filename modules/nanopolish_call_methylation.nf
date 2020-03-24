process call_methylation {
  publishDir "${params.outdir}/${prefix}/methylation", mode: 'copy'
  container = 'fmalmeida/bacannot:latest'
  tag "Methylated sites are being calculated with Nanopolish"

  input:
  tuple val(prefix), file(draft), file(reads), file(fast5)

  output:
  // Grab all outputs
  file "*_calls.tsv" optional true
  file "*_frequency.tsv" optional true
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("cpg_frequency.bedGraph") optional true
  tuple val(prefix), file("gpc_frequency.bedGraph") optional true
  tuple val(prefix), file("dam_frequency.bedGraph") optional true
  tuple val(prefix), file("dcm_frequency.bedGraph") optional true
  tuple val(prefix), file("chr.sizes") optional true

  when:
  // When an entry does not exist, it is created as 'input'
  if ("${fast5_dir}" != 'input' && "${reads}" != 'input.2')

  script:
  fast5_dir = fast5.baseName
  """
  # Index Our Fast5 Data
  nanopolish index -d ${fast5_dir} ${reads} ;

  # Map Our Indexed Reads to Our Genome
  minimap2 -a -x map-ont ${draft} ${reads} | samtools sort -T tmp -o reads_output.sorted.bam ;
  samtools index reads_output.sorted.bam ;

  # Call Methylation
  ## cpg , gpc , dam and dcm
  nanopolish call-methylation -r ${reads} -b reads_output.sorted.bam -g ${draft} --methylation cpg > cpg_calls.tsv ;
	nanopolish call-methylation -r ${reads} -b reads_output.sorted.bam -g ${draft} --methylation gpc > gpc_calls.tsv ;
	nanopolish call-methylation -r ${reads} -b reads_output.sorted.bam -g ${draft} --methylation dam > dam_calls.tsv ;
	nanopolish call-methylation -r ${reads} -b reads_output.sorted.bam -g ${draft} --methylation dcm > dcm_calls.tsv ;

  # Calculate Methylation Frequencies
  ## cpg , gpc , dam and dcm
  /work/nanopolish_scripts/calculate_methylation_frequency.py cpg_calls.tsv > cpg_frequency.tsv ;
  /work/nanopolish_scripts/calculate_methylation_frequency.py gpc_calls.tsv > gpc_frequency.tsv ;
  /work/nanopolish_scripts/calculate_methylation_frequency.py dam_calls.tsv > dam_frequency.tsv ;
  /work/nanopolish_scripts/calculate_methylation_frequency.py dcm_calls.tsv > dcm_frequency.tsv ;

  # Transform These TSV files into bedGraph
  ## cpg
  [ ! -s cpg_frequency.tsv ] || grep -v "start" cpg_frequency.tsv | \
  awk '{ print \$1 "\t" \$2 "\t" \$3 "\t" \$7 }' > cpg_frequency.bedGraph ;
  ## gpc
  [ ! -s gpc_frequency.tsv ] || grep -v "start" gpc_frequency.tsv | \
  awk '{ print \$1 "\t" \$2 "\t" \$3 "\t" \$7 }' > gpc_frequency.bedGraph ;
  ## dam
  [ ! -s dam_frequency.tsv ] || grep -v "start" dam_frequency.tsv |
  awk '{ print \$1 "\t" \$2 "\t" \$3 "\t" \$7 }' > dam_frequency.bedGraph ;
  ## dcm
  [ ! -s dcm_frequency.tsv ] || grep -v "start" dcm_frequency.tsv |
  awk '{ print \$1 "\t" \$2 "\t" \$3 "\t" \$7 }' > dcm_frequency.bedGraph ;

  # Create Contig Sizes File
  seqtk comp ${draft} | awk '{ print \$1 "\t" \$2 }' > chr.sizes
  """
  //else
  //"""
  //echo "Process was skipped"
  //"""
}
