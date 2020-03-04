process compute_gc {
  container 'fmalmeida/bacannot:latest'

  input:
  file 'input.fasta'

  output:
  file "input_GC_500_bps.sorted.bedGraph" // file containing gc values
  file "input.sizes" into gc_sizes_jbrowse // file containing chr sizes

  """
  # Index
  samtools faidx input.fasta ;
  # Take Sizes
  cut -f 1,2 input.fasta.fai > input.sizes ;
  # Create sliding window
  bedtools makewindows -g input.sizes -w 500 > input_500_bps.bed ;
  # Compute GC
  bedtools nuc -fi input.fasta -bed input_500_bps.bed > input_500_bps_nuc.txt ;
  # Create bedGraph
  awk 'BEGIN{FS="\\t"; OFS="\\t"} FNR > 1 { print \$1,\$2,\$3,\$5 }' input_500_bps_nuc.txt > input_GC_500_bps.bedGraph
  # Sort
  bedtools sort -i input_GC_500_bps.bedGraph > input_GC_500_bps.sorted.bedGraph
  """
}
