process call_methylation {
  if (params.fast5_dir && params.fastq_reads) {
  publishDir "${outDir}/methylation", mode: 'copy' }
  container 'fmalmeida/bacannot:latest'
  x = ( params.fast5_dir && params.fastq_reads
        ? "Methylated sites are being calculated"
        : "Process was skipped by the user")
  tag { x }

  input:
  file 'input.fa' from renamed_genome
  file 'reads.fq' from nanopolish_lreads
  file fast5
  //val fast5_dir from fast5_dir

  output:
  file "*_calls.tsv" optional true
  file "*_frequency.tsv" optional true
  file "cpg_frequency.bedGraph" into cpg_bedGraph
  file "gpc_frequency.bedGraph" into gpc_bedGraph
  file "dam_frequency.bedGraph" into dam_bedGraph
  file "dcm_frequency.bedGraph" into dcm_bedGraph
  file "chr.sizes"  into chr_sizes

  script:
  fast5_dir = fast5.baseName
  if (params.fast5_dir && params.fastq_reads)
  """
  # Index Our Fast5 Data
  nanopolish index -d "${fast5_dir}" reads.fq ;
  # Map Our Indexed Reads to Our Genome
  minimap2 -a -x map-ont input.fa reads.fq | samtools sort -T tmp -o reads_output.sorted.bam ;
  samtools index reads_output.sorted.bam ;
  # Call Methylation
  ## cpg , gpc , dam and dcm
  nanopolish call-methylation -r reads.fq -b reads_output.sorted.bam -g input.fa --methylation cpg > cpg_calls.tsv ;
	nanopolish call-methylation -r reads.fq -b reads_output.sorted.bam -g input.fa --methylation gpc > gpc_calls.tsv ;
	nanopolish call-methylation -r reads.fq -b reads_output.sorted.bam -g input.fa --methylation dam > dam_calls.tsv ;
	nanopolish call-methylation -r reads.fq -b reads_output.sorted.bam -g input.fa --methylation dcm > dcm_calls.tsv ;
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
  seqtk comp input.fa | awk '{ print \$1 "\t" \$2 }' > chr.sizes
  """
  else
  """
  touch cpg_frequency.bedGraph gpc_frequency.bedGraph dam_frequency.bedGraph dcm_frequency.bedGraph chr.sizes
  """
}
