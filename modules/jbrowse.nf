process jbrowse {
  publishDir "${params.outdir}/${prefix}/jbrowse/", mode: 'copy'
  label 'jbrowse'
  tag "Creating the genome browser with JBrowse"

  input:
  tuple val(prefix), file(gff), file(draft), file("prokka_gff"), file(mlst), file(barrnap),
        file(gc_bedGraph), file(gc_chrSizes), file(kofamscan), file(vfdb),
        file(victors), file(amrfinder), file(rgi), file(iceberg), file(phast),
        file(phigaro), file(genomic_islands), file("methylation"), file("chr.sizes")

  output:
  file "*"

  """
  # Get JBrowse Files in working directory
  cp -R /work/jbrowse/* . ;

  # Form -fat FASTA file for JBROWSE
  prepare-refseqs.pl --fasta $draft --key \"${prefix}\" --out \"data\" ;

  # Add GC content Track
  bedGraphToBigWig $gc_bedGraph $gc_chrSizes data/GC_content.bw ;
  add-bw-track.pl --bw_url GC_content.bw --plot --label \"GC Content\" --key \"GC Content\" \
  --category \"GC Content\" --pos_color darkgray ;

  # Add track with all features
  flatfile-to-json.pl --gff $gff --key \"All features\" --trackType CanvasFeatures \
  --trackLabel \"${prefix} annotated features\" --out \"data\" --nameAttributes \"Name,ID,product\";

  # Add tRNA track
  [ \$(grep "tRNA" ${gff} | wc -l) -eq 0 ] || grep "tRNA" ${gff} > tRNAs.gff ;
  [ \$(grep "tRNA" ${gff} | wc -l) -eq 0 ] || flatfile-to-json.pl --gff tRNAs.gff --key \"tRNA Sequences\" --nameAttributes \"Name,ID,product\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} tRNA sequences\" \
  --config '{ "style": { "color": "darkgreen" }, "displayMode": "compact" }' --out \"data\" ;

  remove-track.pl --trackLabel \"${prefix} tRNA sequences\" --dir data &> /tmp/error
  [ \$(grep "tRNA" ${gff} | wc -l) -eq 0 ] || echo \' { \"compress\" : 0, \
                                 \"displayMode\" : \"compact\", \
                                 \"key\" : \"tRNA Sequences\", \
                                 \"label\" : \"${prefix} tRNA sequences\", \
                                 \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                 \"style\" : { \"className\" : \"feature\", \"color\": \"darkgreen\" }, \
                                 \"trackType\" : \"CanvasFeatures\", \
                                 \"type\" : \"CanvasFeatures\", \
                                 \"nameAttributes\" : \"Name,ID,product\",
                                 \"urlTemplate\" : \"tracks/${prefix} tRNA sequences/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json
  [ \$(grep "tRNA" ${gff} | wc -l) -eq 0 ] || rm -f tRNAs.gff ;

  # Add track without hypothetical features|proteins
  [ \$(grep -v "hypothetical" ${gff} | wc -l) -eq 0 ] || grep -v "hypothetical" ${gff} > no_hypothetical ;
  [ \$(grep -v "hypothetical" ${gff} | wc -l) -eq 0 ] || flatfile-to-json.pl --gff no_hypothetical --key \"Not hypothetical features\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} not hypothetical features\" --out \"data\" \
  --nameAttributes \"Name,ID,product\" && rm -f no_hypothetical ;

  # Add track with all hypothetical features|proteins
  [ \$(grep "hypothetical" ${gff} | wc -l) -eq 0 ] || grep "hypothetical" ${gff} > hypothetical ;
  [ \$(grep "hypothetical" ${gff} | wc -l) -eq 0 ] || flatfile-to-json.pl --gff hypothetical --key \"Only hypothetical features\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} only hypothetical features\" \
  --out \"data\" --nameAttributes \"Name,ID,product\" && rm -f hypothetical ;

  # Add track with all transposases
  [ \$(grep "transposase" ${gff} | wc -l) -eq 0 ] || grep "transposase" ${gff} > transposase ;
  [ \$(grep "transposase" ${gff} | wc -l) -eq 0 ] || flatfile-to-json.pl --gff transposase --key \"Transposases\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} only transposases\" --out \"data\" --nameAttributes \"Name,ID,product\" && rm -f transposase ;

  # Add track with virulence features
  [ \$(grep -e "Virulence" -e "virulence" ${gff} | wc -l) -eq 0 ] || grep -e "Virulence" -e "virulence" ${gff} > virulence ;
  [ \$(grep -e "Virulence" -e "virulence" ${gff} | wc -l) -eq 0 ] || flatfile-to-json.pl --gff virulence --key \"Virulence features from VFDB and Victors databases\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} virulence features from VFDB and Victors databases\" \
  --out \"data\" --nameAttributes \"VFDB_Target,VFDB_Product,Victors_Target,Victors_Product,Name,ID,product\" && rm -f virulence ;

  # Add track with resistance features
  [ \$(grep -e "Resistance" -e "resistance" ${gff} | wc -l) -eq 0 ] || grep -e "Resistance" -e "resistance" ${gff} > resistance ;
  [ \$(grep -e "Resistance" -e "resistance" ${gff} | wc -l) -eq 0 ] || flatfile-to-json.pl --gff resistance --key \"Resistance features from any source\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} resistance features from all sources\" \
  --out \"data\" --nameAttributes \"Name,ID,product\" && rm -f resistance ;

  # Add track with resistance AMRFinder features
  [ \$(grep "AMRFinderPlus" ${gff} | wc -l) -eq 0 ] || grep "AMRFinderPlus" ${gff} > amrfinder ;
  [ \$(grep "AMRFinderPlus" ${gff} | wc -l) -eq 0 ] || flatfile-to-json.pl --gff amrfinder --key \"AMRFinderPLus AMR features\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} resistance features from AMRFinderPlus\" \
  --out \"data\" --nameAttributes \"Gene_Name,Gene_Product,Name,ID,product\" && rm -f amrfinder ;

  # Add track with resistance RGI features
  [ \$(grep "CARD" ${gff} | wc -l) -eq 0 ] || grep "CARD" ${gff} > rgi ;
  [ \$(grep "CARD" ${gff} | wc -l) -eq 0 ] || flatfile-to-json.pl --gff rgi --key \"CARD-RGI AMR features\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} resistance features from CARD-RGI\" \
  --out \"data\" --nameAttributes \"CARD_Target,CARD_Resistance_Mechanism,Name,ID,product\" && rm -f rgi ;

  # Add track with ICEs features
  [ \$(grep "ICEberg" ${gff} | wc -l) -eq 0 ] || grep "ICEberg" ${gff} > ices ;
  [ \$(grep "ICEberg" ${gff} | wc -l) -eq 0 ] || flatfile-to-json.pl --gff ices --key \"ICE genes\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} genes of integrative and conjugative elements from ICEberg database\" --out \"data\" \
  --nameAttributes \"ICEberg_Target,ICEberg_Product,Name,ID,product\" && rm -f ices ;

  # Add track with prophage features
  [ \$(grep "PHAST" ${gff} | wc -l) -eq 0 ] || grep "PHAST" ${gff} > prophage ;
  [ \$(grep "PHAST" ${gff} | wc -l) -eq 0 ] || flatfile-to-json.pl --gff prophage --key \"Prophage genes\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} prophage genes from PHAST database\" --out \"data\" \
  --nameAttributes \"PHAST_Target,PHAST_Product,Name,ID,product\" && rm -f prophage ;

  # Add track with prophage sequences
  [ -s ${phigaro} ] || flatfile-to-json.pl --bed ${phigaro} --key \"Putative prophage Sequences\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} putative prophage sequences predicted by Phigaro\" \
  --config '{ "style": { "color": "blue" }, "displayMode": "compact" }' --out \"data\" ;
  remove-track.pl --trackLabel \"${prefix} putative prophage sequences predicted by Phigaro\" --dir data &> /tmp/error

  [ -s ${phigaro} ] || echo \' { \"compress\" : 0, \
                                     \"displayMode\" : \"compact\", \
                                     \"key\" : \"Putative prophage Sequences\", \
                                     \"label\" : \"${prefix} putative prophage sequences predicted by Phigaro\", \
                                     \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                     \"style\" : { \"className\" : \"feature\", \"color\": \"blue\" }, \
                                     \"trackType\" : \"CanvasFeatures\", \
                                     \"type\" : \"CanvasFeatures\", \
                                     \"urlTemplate\" : \"tracks/${prefix} putative prophage sequences predicted by Phigaro/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json

  # Add track with GIs
  [ -s ${genomic_islands} ] || flatfile-to-json.pl --bed ${genomic_islands} --key \"Genomic Islands\" \
                              --trackType CanvasFeatures --trackLabel \"${prefix} genomic islands\" \
                              --config '{ "style": { "color": "cyan" }, "displayMode": "compact" }' --out \"data\" ;

  # Remove track for configuration
  [ -s ${genomic_islands} ] || remove-track.pl --trackLabel \"${prefix} genomic islands\" --dir data &> /tmp/error
  # Re-create
  [ -s ${genomic_islands} ] || echo \' { \"compress\" : 0, \
                                    \"displayMode\" : \"compact\", \
                                    \"key\" : \"Genomic Islands\", \
                                    \"label\" : \"${prefix} genomic islands\", \
                                    \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                    \"style\" : { \"className\" : \"feature\", \"color\": \"cyan\" }, \
                                    \"trackType\" : \"CanvasFeatures\", \
                                    \"type\" : \"CanvasFeatures\", \
                                    \"urlTemplate\" : \"tracks/${prefix} genomic islands/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json

  # Add track with rRNA sequences
  [ -s ${barrnap} ] || flatfile-to-json.pl --gff ${barrnap} --key \"rRNA Sequences\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} rRNA sequences\" \
  --config '{ "style": { "color": "blue" }, "displayMode": "compact" }' --out \"data\" ;

  remove-track.pl --trackLabel \"${prefix} rRNA sequences\" --dir data &> /tmp/error
  [ -s ${barrnap} ] || echo \' { \"compress\" : 0, \
                                 \"displayMode\" : \"compact\", \
                                 \"key\" : \"rRNA Sequences\", \
                                 \"label\" : \"${prefix} rRNA sequences\", \
                                 \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                 \"style\" : { \"className\" : \"feature\", \"color\": \"blue\" }, \
                                 \"trackType\" : \"CanvasFeatures\", \
                                 \"type\" : \"CanvasFeatures\", \
                                 \"urlTemplate\" : \"tracks/${prefix} rRNA sequences/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json

  # Form -fat bedGraphs
  ## cpg
  [ ! -s methylation ] || bedGraphToBigWig methylation chr.sizes data/methylation.bw ;
  [ ! -s data/methylation.bw ] || add-bw-track.pl --bw_url methylation.bw --plot --label "5mC (CpG) Methylations" --key "5mC (CpG) Methylations" --category "Methylations" \
  --pos_color blue ;
  """
}
