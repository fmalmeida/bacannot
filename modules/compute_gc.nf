process compute_gc {
  container = 'fmalmeida/bacannot:latest'
  tag "Calculating genome GC with bedtools"

  input:
  tuple val(prefix), file(genome)

  output:
  // Outputs must be linked to each prefix (tag)
  tuple val(prefix), file("input_GC_500_bps.sorted.bedGraph"), file("input.sizes")

  """
  # Index
  samtools faidx $genome ;
  # Take Sizes
  cut -f 1,2 ${genome}.fai > input.sizes ;
  # Create sliding window
  bedtools makewindows -g input.sizes -w 500 > input_500_bps.bed ;
  # Compute GC
  bedtools nuc -fi ${genome} -bed input_500_bps.bed > input_500_bps_nuc.txt ;
  # Create bedGraph
  awk 'BEGIN{FS="\\t"; OFS="\\t"} FNR > 1 { print \$1,\$2,\$3,\$5 }' input_500_bps_nuc.txt > input_GC_500_bps.bedGraph
  # Sort
  bedtools sort -i input_GC_500_bps.bedGraph > input_GC_500_bps.sorted.bedGraph
  """
}
