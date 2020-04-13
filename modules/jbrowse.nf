process jbrowse {
  publishDir "${params.outdir}/${prefix}/jbrowse/", mode: 'copy'
  container 'fmalmeida/bacannot:jbrowse'

  input:
  tuple val(prefix), file(gff), file(draft), file("prokka_gff"), file(mlst), file(barrnap),
        file(gc_bedGraph), file(gc_chrSizes), file(kofamscan), file(vfdb), file(victors),
        file(amrfinder), file(rgi), file(iceberg), file(phast), file(phigaro),
        file(genomic_islands), file("cpg"), file("gpc"), file("dam"), file("dcm"), file("chr.sizes")

  output:
  file "*"

  """
  # Get JBrowse Files in working directory
  cp -R /work/jbrowse/* . ;

  # Format FASTA file for JBROWSE
  prepare-refseqs.pl --fasta $draft --key \"${prefix}\" --out \"data\" ;

  # Add GC content Track
  bedGraphToBigWig $gc_bedGraph $gc_chrSizes data/GC_content.bw ;
  add-bw-track.pl --bw_url GC_content.bw --plot --label \"GC Content\" --key \"GC Content\" \
  --category \"GC Content\" --pos_color darkgray ;

  # Add track with all features
  flatfile-to-json.pl --gff $gff --key \"All features\" --trackType CanvasFeatures \
  --trackLabel \"${prefix} annotated features\" --out \"data\" --nameAttributes \"Name,ID,product\";

  # Add tRNA track
  awk '{ if (\$3 == "tRNA" ) print }' ${gff} > tRNAs.gff ;
  [ -s tRNAs.gff ] || flatfile-to-json.pl --gff tRNAs.gff --key \"tRNA Sequences\" --nameAttributes \"Name,ID,product\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} tRNA sequences\" \
  --config '{ "style": { "color": "darkgreen" }, "displayMode": "compact" }' --out \"data\" ;

  remove-track.pl --trackLabel \"${prefix} tRNA sequences\" --dir data &> /tmp/error
  [ -s tRNAs.gff ] || echo \' { \"compress\" : 0, \
                                 \"displayMode\" : \"compact\", \
                                 \"key\" : \"tRNA Sequences\", \
                                 \"label\" : \"${prefix} tRNA sequences\", \
                                 \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                 \"style\" : { \"className\" : \"feature\", \"color\": \"darkgreen\" }, \
                                 \"trackType\" : \"CanvasFeatures\", \
                                 \"type\" : \"CanvasFeatures\", \
                                 \"nameAttributes\" : \"Name,ID,product\",
                                 \"urlTemplate\" : \"tracks/${prefix} tRNA sequences/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json
  rm tRNAs.gff ;

  # Add track without hypothetical features|proteins
  grep -v "hypothetical" ${gff} > no_hypothetical ;
  [ -s no_hypothetical ] || flatfile-to-json.pl --gff no_hypothetical --key \"Not hypothetical features\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} not hypothetical features\" --out \"data\" \
  --nameAttributes \"Name,ID,product\" ;
  rm no_hypothetical ;

  # Add track with all hypothetical features|proteins
  grep "hypothetical" ${gff} > hypothetical ;
  [ -s hypothetical ] || flatfile-to-json.pl --gff hypothetical --key \"Only hypothetical features\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} only hypothetical features\" \
  --out \"data\" --nameAttributes \"Name,ID,product\" ;
  rm hypothetical ;

  # Add track with all transposases
  grep "transposase" ${gff} > transposase ;
  [ -s transposase ] || flatfile-to-json.pl --gff transposase --key \"Transposases\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} only transposases\" --out \"data\" --nameAttributes \"Name,ID,product\" ;
  rm transposase ;

  # Add track with virulence features
  grep "Virulence" ${gff} > virulence ;
  [ -s virulence ] || flatfile-to-json.pl --gff virulence --key \"Virulence features\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} virulence features\" \
  --out \"data\" --nameAttributes \"VFDB_Target,Victors_Target,Name,ID,product\" ;
  rm virulence ;

  # Add track with resistance features
  grep "Resistance" ${gff} > resistance ;
  [ -s resistance ] || flatfile-to-json.pl --gff resistance --key \"All resistance features from all sources\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} resistance features from all sources\" \
  --out \"data\" --nameAttributes \"Name,ID,product\" ;
  rm resistance ;

  # Add track with resistance AMRFinder features
  grep "AMRFinderPlus" ${gff} > amrfinder ;
  [ -s amrfinder ] || flatfile-to-json.pl --gff amrfinder --key \"AMRFinderPLus features\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} resistance features from AMRFinderPlus\" \
  --out \"data\" --nameAttributes \"Gene_Name,Gene_Product,Name,ID,product\" ;
  rm amrfinder ;

  # Add track with resistance RGI features
  grep "CARD" ${gff} > rgi ;
  [ -s rgi ] || flatfile-to-json.pl --gff rgi --key \"CARD-RGI features\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} resistance features from CARD-RGI\" \
  --out \"data\" --nameAttributes \"CARD_Target,CARD_Resistance_Mechanism,Name,ID,product\" ;
  rm rgi ;

  # Add track with ICEs features
  grep "ICEberg" ${gff} > ices ;
  [ -s ices ] || flatfile-to-json.pl --gff ices --key \"ICE features\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} integrative and conjugative elements\" --out \"data\" \
  --nameAttributes \"ICEberg_Target,Name,ID,product\" ;
  rm ices ;

  # Add track with prophage features
  grep "PHAST" ${gff} > prophage ;
  [ -s prophage ] || flatfile-to-json.pl --gff prophage --key \"Prophage features\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} prophage features (PHAST database)\" --out \"data\" \
  --nameAttributes \"PHAST_Target,Name,ID,product\" ;
  rm prophage ;

  # Add track with prophage sequences
  [ -s ${phigaro} ] || flatfile-to-json.pl --bed ${phigaro} --key \"Prophage Sequences\" \
  --trackType CanvasFeatures --trackLabel \"${prefix} prophage sequences (Phigaro prediction)\" \
  --config '{ "style": { "color": "blue" }, "displayMode": "compact" }' --out \"data\" ;
  remove-track.pl --trackLabel \"${prefix} prophage sequences (Phigaro prediction)\" --dir data &> /tmp/error

  [ -s ${phigaro} ] || echo \' { \"compress\" : 0, \
                                     \"displayMode\" : \"compact\", \
                                     \"key\" : \"Prophage Sequences\", \
                                     \"label\" : \"${prefix} prophage sequences (Phigaro prediction)\", \
                                     \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                     \"style\" : { \"className\" : \"feature\", \"color\": \"blue\" }, \
                                     \"trackType\" : \"CanvasFeatures\", \
                                     \"type\" : \"CanvasFeatures\", \
                                     \"urlTemplate\" : \"tracks/${prefix} prophage sequences/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json

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

  # Add track with efflux features
  #[ -s efflux ] || flatfile-to-json.pl --gff efflux --key \"Efflux [pumps] features\" \
  #--trackType CanvasFeatures --trackLabel \"${prefix} efflux features\" --out \"data\" ;
  # Add track with conjugation features
  #[ -s conjugation ] || flatfile-to-json.pl --gff conjugation --key \"Conjugation related features\" \
  #--trackType CanvasFeatures --trackLabel \"${prefix} Conjugation related features\" --out \"data\" ;

  # Format bedGraphs
  ## cpg
  [ -s cpg ] || bedGraphToBigWig cpg chr.sizes data/cpg.bw ;
  ## gpc
  [ -s gpc ] || bedGraphToBigWig gpc chr.sizes data/gpc.bw ;
  ## dam
  [ -s dam ] || bedGraphToBigWig dam chr.sizes data/dam.bw ;
  ## dcm
  [ -s dcm ] || bedGraphToBigWig dcm chr.sizes data/dcm.bw ;
  # Add BigWigs
  ## cpg
  [ -s data/cpg.bw ] || add-bw-track.pl --bw_url cpg.bw --plot --label "CpG Methylations" --key "CpG Methylations" --category "Methylations" \
  --pos_color blue ;
  ## gpc
  [ -s data/gpc.bw ] || add-bw-track.pl --bw_url gpc.bw --plot --label "GpC Methylations" --key "GpC Methylations" --category "Methylations" \
  --pos_color purple ;
  ## dam
  [ -s data/dam.bw ] || add-bw-track.pl --bw_url dam.bw --plot --label "Dam Methylations" --key "Dam Methylations" --category "Methylations" \
  --pos_color pink ;
  ## dcm
  [ -s data/dcm.bw ] || add-bw-track.pl --bw_url dcm.bw --plot --label "Dcm Methylations" --key "Dcm Methylations" --category "Methylations" \
  --pos_color cyan ;
  """
}
