process COMPUTE_GC {
  tag "${prefix}"
  label = [ 'misc', 'process_ultralow' ]

  input:
  tuple val(prefix), file(genome)

  output:
  tuple val(prefix), path("input_GC_500_bps.sorted.bedGraph"), path("input.sizes"), emit: results

  script:
  """
  # Index de genome
  samtools faidx $genome ;

  # Get contig sizes
  cut -f 1,2 ${genome}.fai > input.sizes ;

  # Create genome sliding window
  bedtools makewindows -g input.sizes -w 500 > input_500_bps.bed ;

  # Compute GC content
  bedtools nuc -fi ${genome} -bed input_500_bps.bed > input_500_bps_nuc.txt ;

  # Create bedGraph for JBrowse
  awk 'BEGIN{FS="\\t"; OFS="\\t"} FNR > 1 { print \$1,\$2,\$3,\$5 }' input_500_bps_nuc.txt > input_GC_500_bps.bedGraph

  # Sort results
  bedtools sort -i input_GC_500_bps.bedGraph > input_GC_500_bps.sorted.bedGraph
  """
}
