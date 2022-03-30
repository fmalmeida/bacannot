#!/bin/bash

# Bacannot shell script for plotting the Genomic Islands predicted with
# Islandpath-DIMOB. It uses the python package gff-toolbox.
#
# Author: Felipe M. Almeida (almeidafmarques@outlook.com)

# Help
Help()
{

	# Display Help
	echo
	echo "Simple help message for the utilization of this script"
	echo "It takes the Islands predicted and draws them into separated plots"
	echo
	echo "Syntax: draw_gis.sh [-h|i|f|g]"
	echo "options:"
	echo "h		Print this help"
	echo "i		Predicted Islands in BED"
	echo "f		File with the list of GFFs for plot"
	echo "g		Merged GFF for separation"
	echo
}

# Get the options
while getopts "hi:f:g:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
			i) # get bed
				 gis_bed="$OPTARG"
				 ;;
			g) # get GFF
				 in_gff="$OPTARG"
				 ;;
			f) # get input fofn
				 in_fofn="$OPTARG"
				 ;;
   esac
done

# Main

# First step
# Separate the GFF in smaller subsets
grep "Resistance" $in_gff > ./resistance.gff ;
grep "Virulence"  $in_gff > ./virulence.gff ;
grep "Prophage"   $in_gff > ./prophages.gff ;
grep "ICE"        $in_gff > ./ices.gff ;
grep -v -e "Resistance" -e "Virulence" -e "Prophage" -e "ICE" $in_gff > ./remainder.gff

# Second step
# Iterate over Islands BED and plot each one
int=0
while read line ; do
	int=$((int + 1)) ;
	contig=$(echo $line | cut -d " " -f 1) ;
	start=$(echo $line | cut -d " " -f 2) ;
	end=$(echo $line | cut -d " " -f 3) ;
	# ID labelled
	gff-toolbox plot --fofn "$in_fofn" --start "$start" --end "$end" --contig "$contig" -t "${contig}: Genomic Islands ${int}" \
		--feature CDS,rRNA,tRNA --output "plots/id_label/${contig}_GI_${int}.png" --width 25 --height 10

	# Produtc labelled
	gff-toolbox plot --fofn "$in_fofn" --start "$start" --end "$end" --contig "$contig" -t "${contig}: Genomic Islands ${int}" \
		--feature CDS,rRNA,tRNA --output "plots/product_label/${contig}_GI_${int}.png" --identification "product" --width 25 --height 10
done<$gis_bed

# Clean
# rm -f *.gff
