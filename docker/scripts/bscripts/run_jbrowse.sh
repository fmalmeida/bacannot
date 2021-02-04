#!/bin/bash

# Bacannot shell script for plotting the JBrowse genome browser
# in a more reproducible and readable manner
#
# Author: Felipe M. Almeida (almeidafmarques@outlook.com)

# Help
Help()
{

	# Display Help
	echo
	echo "Simple help message for the utilization of this script"
	echo "It takes the jbrowse data path and all the files that shall be plotted from bacannot"
	echo
	echo "Syntax: run_jbrowse.sh [-h|p|g|b|s|f|r|B|P|G|m|S|R]"
	echo "options:"
	echo
	echo "h					Print this help"
	echo "p					Prefix for labelling the features"
	echo "g					Path to genome in FASTA file"
	echo "b					Path to GC content in bedGraph"
	echo "s					Path to file with chr sizes"
	echo "f					Path to complete GFF file (from prokka)"
	echo "r					Path to barrnap GFF"
	echo "B					Path to prophage sequences in BED format from phigaro"
	echo "P					Path to prophage sequences in BED format from phispy"
	echo "G					Path to genomic islands in BED format from IslandPath-DIMOB"
	echo "m					Path to Nanopolish methylation results"
	echo "S					Path to Nanopolish chr sizes"
	echo "R					Path to Resfinder custom GFF"
	echo ""
	echo
}

# Get the options
while getopts "hp:g:b:s:f:r:B:P:G:m:S:R:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
			p) # get genome prefix
				 PREFIX="$OPTARG"
				 ;;
			g) # get genome FASTA
				 GENOME="$OPTARG"
				 ;;
			b) # get GC bedgraph
				 BEDGRAPH="$OPTARG"
				 ;;
			s) # get chr sizes
				 CHRSIZES="$OPTARG"
				 ;;
			f) # get prokka gff
				 PROKKAGFF="$OPTARG"
				 ;;
			r) # get barrnap gff
				 rRNAGFF="$OPTARG"
				 ;;
			B) # get phigaro bed
				 PHIGAROBED="$OPTARG"
				 ;;
			P) # get phispy bed
				 PHISPYBED="$OPTARG"
				 ;;
			G) # get GIs bed
				 GIBED="$OPTARG"
				 ;;
			m) # get nanopolish methylation
				 NANOMETHYL="$OPTARG"
				 ;;
			S) # get nanopolish chr sizes
				 NANOSIZES="$OPTARG"
				 ;;
			R) # get resfinder GFF
				 RESFINDERGFF="$OPTARG"
				 ;;
   esac
done

# Main

# First STEP, index the genome fasta
prepare-refseqs.pl --fasta $GENOME --key "$PREFIX" --out "data" ;

# Add GC content Track
bedGraphToBigWig $BEDGRAPH $CHRSIZES data/GC_content.bw ;
add-bw-track.pl --bw_url GC_content.bw --plot --label "GC Content" --key "GC Content" --category "GC Content" --pos_color darkgray ;

# Add track with all features
flatfile-to-json.pl --gff $PROKKAGFF --key "${PREFIX} all features" --trackType CanvasFeatures \
--trackLabel "${PREFIX} all features" --out "data" --nameAttributes "Name,ID,gene,product";
remove-track.pl --trackLabel "${PREFIX} all features" --dir data &> /tmp/error
echo -E "{ \"compress\" : 0, \
					 \"displayMode\" : \"compact\", \
					 \"key\" : \"${PREFIX} all features\", \
					 \"category\" : \"Generic annotation\", \
					 \"label\" : \"${PREFIX} all features\", \
					 \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
					 \"style\" : { \"className\" : \"feature\", \"color\": \"goldenrod\" }, \
					 \"trackType\" : \"CanvasFeatures\", \
					 \"type\" : \"CanvasFeatures\", \
					 \"nameAttributes\" : \"Name,ID,gene,product\", \
					 \"urlTemplate\" : \"tracks/${PREFIX} all features/{refseq}/trackData.json\" }" | add-track-json.pl  data/trackList.json

# Add tRNA track
[ $(grep "tRNA" $PROKKAGFF  | wc -l) -eq 0 ] || awk '{ if ($3 == "tRNA") print }' $PROKKAGFF  > tRNAs.gff ;
[ $(grep "tRNA" $PROKKAGFF  | wc -l) -eq 0 ] || flatfile-to-json.pl --gff tRNAs.gff --key "${PREFIX} tRNA sequences" --nameAttributes "Name,ID,product" \
--trackType CanvasFeatures --trackLabel "${PREFIX} tRNA sequences" --config '{ "style": { "color": "darkgreen" }, "displayMode": "compact" }' --out "data" ;
remove-track.pl --trackLabel "${PREFIX} tRNA sequences" --dir data &> /tmp/error
[ $(grep "tRNA" $PROKKAGFF  | wc -l) -eq 0 ] || echo -E " { \"compress\" : 0, \
																												 		\"displayMode\" : \"compact\", \
																														\"key\" : \"${PREFIX} tRNA sequences\", \
																														\"category\" : \"Generic annotation\", \
																														\"label\" : \"${PREFIX} tRNA sequences\", \
																														\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																														\"style\" : { \"className\" : \"feature\", \"color\": \"darkgreen\" }, \
																														\"trackType\" : \"CanvasFeatures\", \
																														\"type\" : \"CanvasFeatures\", \
																														\"nameAttributes\" : \"Name,ID,product\", \
																														\"urlTemplate\" : \"tracks/${PREFIX} tRNA sequences/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ $(grep "tRNA" $PROKKAGFF  | wc -l) -eq 0 ] || rm -f tRNAs.gff ;

# Add track with rRNA sequences
[ ! -s $rRNAGFF ] || flatfile-to-json.pl --gff $rRNAGFF --key "${PREFIX} rRNA Sequences" --trackType CanvasFeatures --trackLabel "${PREFIX} rRNA Sequences" \
--config '{ "style": { "color": "blue" }, "displayMode": "compact" }' --out "data" ;
remove-track.pl --trackLabel "${PREFIX} rRNA Sequences" --dir data &> /tmp/error
[ ! -s $rRNAGFF ]  || echo -E " { \"compress\" : 0, \
															 		\"displayMode\" : \"compact\", \
																	\"key\" : \"${PREFIX} rRNA Sequences\", \
																	\"category\" : \"Generic annotation\", \
																	\"label\" : \"${PREFIX} rRNA Sequences\", \
																	\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																	\"style\" : { \"className\" : \"feature\", \"color\": \"blue\" }, \
																	\"trackType\" : \"CanvasFeatures\", \
																	\"type\" : \"CanvasFeatures\", \
																	\"nameAttributes\" : \"Name,ID,product\", \
																	\"urlTemplate\" : \"tracks/${PREFIX} rRNA Sequences/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json

# Add track without hypothetical features|proteins
[ $(grep -v "hypothetical" $PROKKAGFF  | wc -l) -eq 0 ] || grep -v "hypothetical" $PROKKAGFF  > no_hypothetical ;
[ $(grep -v "hypothetical" $PROKKAGFF  | wc -l) -eq 0 ] || flatfile-to-json.pl --gff no_hypothetical --key "${PREFIX} not hypothetical features" \
--trackType CanvasFeatures --trackLabel "${PREFIX} not hypothetical features" --out "data" --nameAttributes "Name,ID,product" && rm -f no_hypothetical ;
remove-track.pl --trackLabel "${PREFIX} not hypothetical features" --dir data &> /tmp/error
[ $(grep -v "hypothetical" $PROKKAGFF  | wc -l) -eq 0 ] || echo -E " {  \"compress\" : 0, \
																																		 		\"displayMode\" : \"compact\", \
																																				\"key\" : \"${PREFIX} not hypothetical features\", \
																																				\"category\" : \"Generic annotation\", \
																																				\"label\" : \"${PREFIX} not hypothetical features\", \
																																				\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																																				\"style\" : { \"className\" : \"feature\", \"color\": \"goldenrod\" }, \
																																				\"trackType\" : \"CanvasFeatures\", \
																																				\"type\" : \"CanvasFeatures\", \
																																				\"nameAttributes\" : \"Name,ID,product\", \
																																				\"urlTemplate\" : \"tracks/${PREFIX} not hypothetical features/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ $(grep -v "hypothetical" $PROKKAGFF  | wc -l) -eq 0 ] || rm -f no_hypothetical ;

# Add track with all hypothetical features|proteins
[ $(grep "hypothetical" $PROKKAGFF  | wc -l) -eq 0 ] || grep "hypothetical" $PROKKAGFF > hypothetical ;
[ $(grep "hypothetical" $PROKKAGFF  | wc -l) -eq 0 ] || flatfile-to-json.pl --gff hypothetical --key "${PREFIX} hypothetical features" \
--trackType CanvasFeatures --trackLabel "${PREFIX} hypothetical features"   --out "data" --nameAttributes "Name,ID,product" && rm -f hypothetical ;
remove-track.pl --trackLabel "${PREFIX} hypothetical features" --dir data &> /tmp/error
[ $(grep "hypothetical" $PROKKAGFF  | wc -l) -eq 0 ] || echo -E " {  \"compress\" : 0, \
																																		 		\"displayMode\" : \"compact\", \
																																				\"key\" : \"${PREFIX} hypothetical features\", \
																																				\"category\" : \"Generic annotation\", \
																																				\"label\" : \"${PREFIX} hypothetical features\", \
																																				\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																																				\"style\" : { \"className\" : \"feature\", \"color\": \"darkgray\" }, \
																																				\"trackType\" : \"CanvasFeatures\", \
																																				\"type\" : \"CanvasFeatures\", \
																																				\"nameAttributes\" : \"Name,ID,product\", \
																																				\"urlTemplate\" : \"tracks/${PREFIX} hypothetical features/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ $(grep "hypothetical" $PROKKAGFF  | wc -l) -eq 0 ] || rm -f hypothetical ;

# Add track with all transposases
[ $(grep "transposase" $PROKKAGFF  | wc -l) -eq 0 ] || grep "transposase" $PROKKAGFF  > transposase ;
[ $(grep "transposase" $PROKKAGFF  | wc -l) -eq 0 ] || flatfile-to-json.pl --gff transposase --key "${PREFIX} transposases" --trackType CanvasFeatures \
--trackLabel "${PREFIX} transposases" --out "data" --nameAttributes "Name,ID,product" && rm -f transposase ;
remove-track.pl --trackLabel "${PREFIX} transposases" --dir data &> /tmp/error
[ $(grep "transposase" $PROKKAGFF  | wc -l) -eq 0 ] || echo -E " {  				\"compress\" : 0, \
																																		 		\"displayMode\" : \"compact\", \
																																				\"key\" : \"${PREFIX} transposases\", \
																																				\"category\" : \"Generic annotation\", \
																																				\"label\" : \"${PREFIX} transposases\", \
																																				\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																																				\"style\" : { \"className\" : \"feature\", \"color\": \"#e7a134\" }, \
																																				\"trackType\" : \"CanvasFeatures\", \
																																				\"type\" : \"CanvasFeatures\", \
																																				\"nameAttributes\" : \"Name,ID,product\", \
																																				\"urlTemplate\" : \"tracks/${PREFIX} transposases/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ $(grep "transposase" $PROKKAGFF  | wc -l) -eq 0 ] || rm -f transposase ;

# Add track with virulence features
[ $(grep -e "Virulence" -e "virulence" $PROKKAGFF | wc -l) -eq 0 ] || grep -e "Virulence" -e "virulence" $PROKKAGFF  > virulence ;
[ $(grep -e "Virulence" -e "virulence" $PROKKAGFF | wc -l) -eq 0 ] || flatfile-to-json.pl --gff virulence --key "${PREFIX} all virulence features" \
--trackType CanvasFeatures --trackLabel "${PREFIX} all virulence features" --out "data" --nameAttributes "VFDB_Target,VFDB_Product,Victors_Target,Victors_Product,Name,ID,product"
remove-track.pl --trackLabel "${PREFIX} all virulence features" --dir data &> /tmp/error
[ $(grep -e "Virulence" -e "virulence" $PROKKAGFF | wc -l) -eq 0 ] || echo -E " { \"compress\" : 0, \
																																						 		\"displayMode\" : \"compact\", \
																																								\"key\" : \"${PREFIX} all virulence features\", \
																																								\"category\" : \"Virulence annotation\", \
																																								\"label\" : \"${PREFIX} all virulence features\", \
																																								\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																																								\"style\" : { \"className\" : \"feature\", \"color\": \"darkred\" }, \
																																								\"trackType\" : \"CanvasFeatures\", \
																																								\"type\" : \"CanvasFeatures\", \
																																								\"nameAttributes\" : \"VFDB_Target,VFDB_Product,Victors_Target,Victors_Product,Name,ID,product\", \
																																								\"urlTemplate\" : \"tracks/${PREFIX} all virulence features/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ $(grep -e "Virulence" -e "virulence" $PROKKAGFF | wc -l) -eq 0 ] || rm -f virulence ;

## VFDB
[ $(grep -e "VFDB" $PROKKAGFF | wc -l) -eq 0 ] || grep "VFDB" $PROKKAGFF > vfdb ;
[ $(grep -e "VFDB" $PROKKAGFF | wc -l) -eq 0 ] || flatfile-to-json.pl --gff vfdb --key "${PREFIX} VFDB virulence features" \
--trackType CanvasFeatures --trackLabel "${PREFIX} VFDB virulence features" --out "data" --nameAttributes "VFDB_Target,VFDB_Product,Name,ID,product"
remove-track.pl --trackLabel "${PREFIX} VFDB virulence features" --dir data &> /tmp/error
[ $(grep -e "VFDB" $PROKKAGFF | wc -l) -eq 0 ] || echo -E " { \"compress\" : 0, \
																													 		\"displayMode\" : \"compact\", \
																															\"key\" : \"${PREFIX} VFDB virulence features\", \
																															\"category\" : \"Virulence annotation\", \
																															\"label\" : \"${PREFIX} VFDB virulence features\", \
																															\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																															\"style\" : { \"className\" : \"feature\", \"color\": \"darkred\" }, \
																															\"trackType\" : \"CanvasFeatures\", \
																															\"type\" : \"CanvasFeatures\", \
																															\"nameAttributes\" : \"VFDB_Target,VFDB_Product,Name,ID,product\", \
																															\"urlTemplate\" : \"tracks/${PREFIX} VFDB virulence features/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ $(grep -e "VFDB" $PROKKAGFF | wc -l) -eq 0 ] || rm -f vfdb ;

## Victors
[ $(grep -e "Victors" $PROKKAGFF | wc -l) -eq 0 ] || grep "Victors" $PROKKAGFF  > victors ;
[ $(grep -e "Victors" $PROKKAGFF | wc -l) -eq 0 ] || flatfile-to-json.pl --gff victors --key "${PREFIX} Victors virulence features" \
--trackType CanvasFeatures --trackLabel "${PREFIX} Victors virulence features" --out "data" --nameAttributes "Victors_Target,Victors_Product,Name,ID,product"
remove-track.pl --trackLabel "${PREFIX} Victors virulence features" --dir data &> /tmp/error
[ $(grep -e "Victors" $PROKKAGFF | wc -l) -eq 0 ] || echo -E " {  \"compress\" : 0, \
																															 		\"displayMode\" : \"compact\", \
																																	\"key\" : \"${PREFIX} Victors virulence features\", \
																																	\"category\" : \"Virulence annotation\", \
																																	\"label\" : \"${PREFIX} Victors virulence features\", \
																																	\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																																	\"style\" : { \"className\" : \"feature\", \"color\": \"darkred\" }, \
																																	\"trackType\" : \"CanvasFeatures\", \
																																	\"type\" : \"CanvasFeatures\", \
																																	\"nameAttributes\" : \"Victors_Target,Victors_Product,Name,ID,product\", \
																																	\"urlTemplate\" : \"tracks/${PREFIX} Victors virulence features/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ $(grep -e "Victors" $PROKKAGFF | wc -l) -eq 0 ] || rm -f victors ;

# Add track with resistance features
[ $(grep -e "Resistance" -e "resistance" $PROKKAGFF | wc -l) -eq 0 ] || grep -e "Resistance" -e "resistance" $PROKKAGFF> resistance ;
[ $(grep -e "Resistance" -e "resistance" $PROKKAGFF | wc -l) -eq 0 ] || flatfile-to-json.pl --gff resistance --key "${PREFIX} all resistance features" \
--trackType CanvasFeatures --trackLabel "${PREFIX} all resistance features" --out "data" --nameAttributes "Name,ID,product" ;
remove-track.pl --trackLabel "${PREFIX} all resistance features" --dir data &> /tmp/error
[ $(grep -e "Resistance" -e "resistance" $PROKKAGFF | wc -l) -eq 0 ] || echo -E " {  \"compress\" : 0, \
																															 		\"displayMode\" : \"compact\", \
																																	\"key\" : \"${PREFIX} all resistance features\", \
																																	\"category\" : \"Resistance annotation\", \
																																	\"label\" : \"${PREFIX} all resistance features\", \
																																	\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																																	\"style\" : { \"className\" : \"feature\", \"color\": \"purple\" }, \
																																	\"trackType\" : \"CanvasFeatures\", \
																																	\"type\" : \"CanvasFeatures\", \
																																	\"nameAttributes\" : \"Name,ID,product\", \
																																	\"urlTemplate\" : \"tracks/${PREFIX} all resistance features/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ $(grep -e "Resistance" -e "resistance" $PROKKAGFF | wc -l) -eq 0 ] || rm -f resistance ;

## AMRFinderPlus
[ $(grep "AMRFinderPlus" $PROKKAGFF | wc -l) -eq 0 ] || grep "AMRFinderPlus" $PROKKAGFF> amrfinder ;
[ $(grep "AMRFinderPlus" $PROKKAGFF | wc -l) -eq 0 ] || flatfile-to-json.pl --gff amrfinder --key "${PREFIX} AMRFinder resistance features" --trackType CanvasFeatures \
--trackLabel "${PREFIX} AMRFinder resistance features" --out "data" --nameAttributes "NDARO_Gene_Name,NDARO_Gene_Product,Name,ID,product" ;
remove-track.pl --trackLabel "${PREFIX} AMRFinder resistance features" --dir data &> /tmp/error
[ $(grep "AMRFinderPlus" $PROKKAGFF | wc -l) -eq 0 ] || echo -E " {  \"compress\" : 0, \
																															 		\"displayMode\" : \"compact\", \
																																	\"key\" : \"${PREFIX} AMRFinder resistance features\", \
																																	\"category\" : \"Resistance annotation\", \
																																	\"label\" : \"${PREFIX} AMRFinder resistance features\", \
																																	\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																																	\"style\" : { \"className\" : \"feature\", \"color\": \"purple\" }, \
																																	\"trackType\" : \"CanvasFeatures\", \
																																	\"type\" : \"CanvasFeatures\", \
																																	\"nameAttributes\" : \"NDARO_Gene_Name,NDARO_Gene_Product,Name,ID,product\", \
																																	\"urlTemplate\" : \"tracks/${PREFIX} AMRFinder resistance features/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ $(grep "AMRFinderPlus" $PROKKAGFF | wc -l) -eq 0 ] || rm -f amrfinder ;

## CARD-RGI
[ $(grep "CARD" $PROKKAGFF | wc -l) -eq 0 ] || grep "CARD" $PROKKAGFF> rgi ;
[ $(grep "CARD" $PROKKAGFF | wc -l) -eq 0 ] || flatfile-to-json.pl --gff rgi --key "${PREFIX} CARD-RGI resistance features" --trackType CanvasFeatures \
--trackLabel "${PREFIX} CARD-RGI resistance features" --out "data" --nameAttributes "CARD_name,CARD_product,Targeted_drug_class,Name,ID,product" ;
remove-track.pl --trackLabel "${PREFIX} CARD-RGI resistance features" --dir data &> /tmp/error
[ $(grep "CARD" $PROKKAGFF | wc -l) -eq 0 ] || echo -E " {  \"compress\" : 0, \
																											 		\"displayMode\" : \"compact\", \
																													\"key\" : \"${PREFIX} CARD-RGI resistance features\", \
																													\"category\" : \"Resistance annotation\", \
																													\"label\" : \"${PREFIX} CARD-RGI resistance features\", \
																													\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																													\"style\" : { \"className\" : \"feature\", \"color\": \"purple\" }, \
																													\"trackType\" : \"CanvasFeatures\", \
																													\"type\" : \"CanvasFeatures\", \
																													\"nameAttributes\" : \"CARD_name,CARD_product,Targeted_drug_class,Name,ID,product\", \
																													\"urlTemplate\" : \"tracks/${PREFIX} CARD-RGI resistance features/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ $(grep "CARD" $PROKKAGFF | wc -l) -eq 0 ] || rm -f rgi ;

## Resfinder
[ ! -s $RESFINDERGFF ] || flatfile-to-json.pl --gff $RESFINDERGFF --key "${PREFIX} Resfinder resistance features" --trackType CanvasFeatures \
--trackLabel "${PREFIX} Resfinder resistance features" --out "data" --nameAttributes "Resfinder_gene,ID,Resfinder_phenotype" ;
remove-track.pl --trackLabel "${PREFIX} Resfinder resistance features" --dir data &> /tmp/error
[ ! -s $RESFINDERGFF ] || echo -E " { \"compress\" : 0, \
																	 		\"displayMode\" : \"compact\", \
																			\"key\" : \"${PREFIX} Resfinder resistance features\", \
																			\"category\" : \"Resistance annotation\", \
																			\"label\" : \"${PREFIX} Resfinder resistance features\", \
																			\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																			\"style\" : { \"className\" : \"feature\", \"color\": \"purple\" }, \
																			\"trackType\" : \"CanvasFeatures\", \
																			\"type\" : \"CanvasFeatures\", \
																			\"nameAttributes\" : \"Resfinder_gene,ID,Resfinder_phenotype\", \
																			\"urlTemplate\" : \"tracks/${PREFIX} Resfinder resistance features/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json

# Add mobile genetic elements
## ICEs
[ $(grep "ICEberg" $PROKKAGFF | wc -l) -eq 0 ] || grep "ICEberg" $PROKKAGFF > ices ;
[ $(grep "ICEberg" $PROKKAGFF | wc -l) -eq 0 ] || flatfile-to-json.pl --gff ices --key "${PREFIX} ICE genes from ICEberg database" --trackType CanvasFeatures \
--trackLabel "${PREFIX} ICE genes from ICEberg database" --out "data" --nameAttributes "ICEberg_Target,ICEberg_Product,Name,ID,product" ;
remove-track.pl --trackLabel "${PREFIX} ICE genes from ICEberg database" --dir data &> /tmp/error
[ $(grep "ICEberg" $PROKKAGFF | wc -l) -eq 0 ] || echo -E " {  \"compress\" : 0, \
																											 		\"displayMode\" : \"compact\", \
																													\"key\" : \"${PREFIX} ICE genes from ICEberg database\", \
																													\"category\" : \"MGEs annotation\", \
																													\"label\" : \"${PREFIX} ICE genes from ICEberg database\", \
																													\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																													\"style\" : { \"className\" : \"feature\", \"color\": \"#6db6d9\" }, \
																													\"trackType\" : \"CanvasFeatures\", \
																													\"type\" : \"CanvasFeatures\", \
																													\"nameAttributes\" : \"ICEberg_Target,ICEberg_Product,Name,ID,product\", \
																													\"urlTemplate\" : \"tracks/${PREFIX} ICE genes from ICEberg database/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ $(grep "ICEberg" $PROKKAGFF | wc -l) -eq 0 ] || rm -f iceberg ;

## PROPHAGES
### PHAST
[ $(grep "PHAST" $PROKKAGFF | wc -l) -eq 0 ] || grep "PHAST" $PROKKAGFF > prophage ;
[ $(grep "PHAST" $PROKKAGFF | wc -l) -eq 0 ] || flatfile-to-json.pl --gff prophage --key "${PREFIX} prophage genes from PHAST database" --trackType CanvasFeatures \
--trackLabel "${PREFIX} prophage genes from PHAST database" --out "data" --nameAttributes "PHAST_Target,PHAST_Product,Name,ID,product" ;
remove-track.pl --trackLabel "${PREFIX} prophage genes from PHAST database" --dir data &> /tmp/error
[ $(grep "PHAST" $PROKKAGFF | wc -l) -eq 0 ] || echo -E " {  \"compress\" : 0, \
																											 		\"displayMode\" : \"compact\", \
																													\"key\" : \"${PREFIX} prophage genes from PHAST database\", \
																													\"category\" : \"MGEs annotation\", \
																													\"label\" : \"${PREFIX} prophage genes from PHAST database\", \
																													\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																													\"style\" : { \"className\" : \"feature\", \"color\": \"#1eacb0\" }, \
																													\"trackType\" : \"CanvasFeatures\", \
																													\"type\" : \"CanvasFeatures\", \
																													\"nameAttributes\" : \"PHAST_Target,PHAST_Product,Name,ID,product\", \
																													\"urlTemplate\" : \"tracks/${PREFIX} prophage genes from PHAST database/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ $(grep "PHAST" $PROKKAGFF | wc -l) -eq 0 ] || rm -f prophage ;

### PHIGARO
[ ! -s $PHIGAROBED ] || flatfile-to-json.pl --bed $PHIGAROBED --key "${PREFIX} putative prophages predicted by phigaro" --trackType CanvasFeatures \
--trackLabel "${PREFIX} putative prophages predicted by phigaro" --config '{ "style": { "color": "#00ffff" }, "displayMode": "compact" }' --out "data" ;
remove-track.pl --trackLabel "${PREFIX} putative prophages predicted by phigaro" --dir data &> /tmp/error
[ ! -s $PHIGAROBED ] || echo -E " {  \"compress\" : 0, \
																											 		\"displayMode\" : \"compact\", \
																													\"key\" : \"${PREFIX} putative prophages predicted by phigaro\", \
																													\"category\" : \"MGEs annotation\", \
																													\"label\" : \"${PREFIX} putative prophages predicted by phigaro\", \
																													\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																													\"style\" : { \"className\" : \"feature\", \"color\": \"#00ffff\" }, \
																													\"trackType\" : \"CanvasFeatures\", \
																													\"type\" : \"CanvasFeatures\", \
																													\"urlTemplate\" : \"tracks/${PREFIX} putative prophages predicted by phigaro/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json

### PHISPY
[ ! -s $PHISPYBED ] || tail -n +2 $PHISPYBED | cut -f 2,3,4 > phispy ;
[ ! -s $PHISPYBED ] || flatfile-to-json.pl --bed phispy --key "${PREFIX} putative prophages predicted by phispy" --trackType CanvasFeatures \
--trackLabel "${PREFIX} putative prophages predicted by phispy" --config '{ "style": { "color": "#1eacb0" }, "displayMode": "compact" }' --out "data" ;
remove-track.pl --trackLabel "${PREFIX} putative prophages predicted by phispy" --dir data &> /tmp/error
[ ! -s $PHISPYBED ] || echo -E " {  \"compress\" : 0, \
																											 		\"displayMode\" : \"compact\", \
																													\"key\" : \"${PREFIX} putative prophages predicted by phispy\", \
																													\"category\" : \"MGEs annotation\", \
																													\"label\" : \"${PREFIX} putative prophages predicted by phispy\", \
																													\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																													\"style\" : { \"className\" : \"feature\", \"color\": \"#1eacb0\" }, \
																													\"trackType\" : \"CanvasFeatures\", \
																													\"type\" : \"CanvasFeatures\", \
																													\"urlTemplate\" : \"tracks/${PREFIX} putative prophages predicted by phispy/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json
[ ! -s $PHISPYBED ] || rm -f phispy ;

## Genomic Islands
[ ! -s $GIBED ] || flatfile-to-json.pl --bed $GIBED --key "${PREFIX} putative genomic islands" --trackType CanvasFeatures \
--trackLabel "${PREFIX} putative genomic islands" --config '{ "style": { "color": "#199db0" }, "displayMode": "compact" }' --out "data" ;
[ ! -s $GIBED ] || remove-track.pl --trackLabel "${PREFIX} putative genomic islands" --dir data &> /tmp/error
[ ! -s $GIBED ] || echo -E " {  \"compress\" : 0, \
														 		\"displayMode\" : \"compact\", \
																\"key\" : \"${PREFIX} putative genomic islands\", \
																\"category\" : \"MGEs annotation\", \
																\"label\" : \"${PREFIX} putative genomic islands\", \
																\"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
																\"style\" : { \"className\" : \"feature\", \"color\": \"#199db0\" }, \
																\"trackType\" : \"CanvasFeatures\", \
																\"type\" : \"CanvasFeatures\", \
																\"urlTemplate\" : \"tracks/${PREFIX} putative genomic islands/{refseq}/trackData.json\" } " | add-track-json.pl  data/trackList.json

# Form -fat bedGraphs
## cpg
[ ! -s $NANOMETHYL ] || bedGraphToBigWig $NANOMETHYL $NANOSIZES data/methylation.bw ;
[ ! -s data/methylation.bw ] || add-bw-track.pl --bw_url methylation.bw --plot --label "5mC (CpG) Methylations" --key "5mC (CpG) Methylations" --category "Methylations" --pos_color "#0e469a" ;

# Finally fix categories order
echo >> jbrowse.conf
echo "[trackSelector]" >> jbrowse.conf
echo "categoryOrder=Reference sequence,Generic annotation,Virulence annotation,Resistance annotation,MGEs annotation,GC Content,Methylations" >> jbrowse.conf
