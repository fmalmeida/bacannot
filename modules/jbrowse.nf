process jbrowse {
  publishDir "${params.outdir}/${prefix}/jbrowse/", mode: 'copy'
  container 'fmalmeida/bacannot:jbrowse'

  input:
  tuple val(prefix), file(draft)
  file gff from final_gff
  file 'GC_content.bedGraph' from gc_content_jbrowse
  file 'GC_content.sizes' from gc_sizes_jbrowse
  file 'rrna.gff' from rrna_gff
  file 'resistance' from resistance_gff
  file 'amrfinder' from amrfinderplus_gff
  file 'rgi' from rgi_gff
  file 'virulence' from virulence_gff
  file 'prophage' from prophage_gff
  file 'prophages.bed' from phigaro_bed
  file 'all_GIs.bed' from predicted_GIs
  file 'ices' from ices_gff
  file 'conjugation' from conjugation_gff
  file 'efflux' from efflux_gff
  file 'no_hypothetical' from noHypothetical_gff
  file 'hypothetical' from hypothetical_gff
  file 'transposase' from transposase_gff
  file 'cpg' from cpg_bedGraph
  file 'gpc' from gpc_bedGraph
  file 'dam' from dam_bedGraph
  file 'dcm' from dcm_bedGraph
  file 'chr.sizes' from chr_sizes

  output:
  file "*" optional true

  """
  # Get Files
  cp -R /work/jbrowse/* . ;
  # Format FASTA file for JBROWSE
  prepare-refseqs.pl --fasta $input --key \"${params.prefix}\" --out \"data\" ;
  # Add GC content Track
  bedGraphToBigWig GC_content.bedGraph GC_content.sizes data/GC_content.bw ;
  add-bw-track.pl --bw_url GC_content.bw --plot --label "GC Content" --key "GC Content" \
  --category "GC Content" --pos_color darkgray ;
  # Add track with all features
  flatfile-to-json.pl --gff $gff --key \"All features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} annotated features\" --out \"data\" ;
  # Add tRNA track
  awk '{ if (\$3 == "tRNA" ) print }' ${gff} > tRNAs.gff ;
  [ ! -s tRNAs.gff ] || flatfile-to-json.pl --gff tRNAs.gff --key \"tRNA Sequences\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} tRNA sequences\" \
  --config '{ "style": { "color": "darkgreen" }, "displayMode": "compact" }' --out \"data\" ;
  remove-track.pl --trackLabel \"${params.prefix} tRNA sequences\" --dir data &> /tmp/error
  [ ! -s tRNAs.gff ] || echo \' { \"compress\" : 0, \
                                 \"displayMode\" : \"compact\", \
                                 \"key\" : \"tRNA Sequences\", \
                                 \"label\" : \"${params.prefix} tRNA sequences\", \
                                 \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                 \"style\" : { \"className\" : \"feature\", \"color\": \"darkgreen\" }, \
                                 \"trackType\" : \"CanvasFeatures\", \
                                 \"type\" : \"CanvasFeatures\", \
                                 \"urlTemplate\" : \"tracks/${params.prefix} tRNA sequences/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json
  # Add track without hypothetical features|proteins
  [ ! -s no_hypothetical ] || flatfile-to-json.pl --gff no_hypothetical --key \"Not hypothetical features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} not hypothetical features\" --out \"data\" ;
  # Add track with all hypothetical features|proteins
  [ ! -s hypothetical ] || flatfile-to-json.pl --gff hypothetical --key \"Only hypothetical features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} only hypothetical features\" --out \"data\" ;
  # Add track with all transposases
  [ ! -s transposase ] || flatfile-to-json.pl --gff transposase --key \"Transposases\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} only transposases\" --out \"data\" ;
  # Add track with virulence features
  [ ! -s virulence ] || flatfile-to-json.pl --gff virulence --key \"Virulence features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} virulence features\" --out \"data\" ;
  # Add track with resistance features
  [ ! -s resistance ] || flatfile-to-json.pl --gff resistance --key \"All resistance features from all sources\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} resistance features from all sources\" --out \"data\" ;
  # Add track with resistance AMRFinder features
  [ ! -s amrfinder ] || flatfile-to-json.pl --gff amrfinder --key \"AMRFinderPLus features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} resistance features from AMRFinderPlus\" --out \"data\" ;
  # Add track with resistance RGI features
  [ ! -s rgi ] || flatfile-to-json.pl --gff rgi --key \"CARD-RGI features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} resistance features from CARD-RGI\" --out \"data\" ;
  # Add track with ICEs features
  [ ! -s ices ] || flatfile-to-json.pl --gff ices --key \"ICE features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} integrative and conjugative elements\" --out \"data\" ;
  # Add track with prophage features
  [ ! -s prophage ] || flatfile-to-json.pl --gff prophage --key \"Prophage features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} prophage features\" --out \"data\" ;
  # Add track with prophage sequences
  [ ! -s prophages.bed ] || flatfile-to-json.pl --bed prophages.bed --key \"Prophage Sequences\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} prophage sequences\" \
  --config '{ "style": { "color": "blue" }, "displayMode": "compact" }' --out \"data\" ;
  remove-track.pl --trackLabel \"${params.prefix} prophage sequences\" --dir data &> /tmp/error
  [ ! -s prophages.bed ] || echo \' { \"compress\" : 0, \
                                     \"displayMode\" : \"compact\", \
                                     \"key\" : \"Prophage Sequences\", \
                                     \"label\" : \"${params.prefix} prophage sequences\", \
                                     \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                     \"style\" : { \"className\" : \"feature\", \"color\": \"blue\" }, \
                                     \"trackType\" : \"CanvasFeatures\", \
                                     \"type\" : \"CanvasFeatures\", \
                                     \"urlTemplate\" : \"tracks/${params.prefix} prophage sequences/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json
  # Add track with GIs
  [ ! -s all_GIs.bed ] || flatfile-to-json.pl --bed all_GIs.bed --key \"Genomic Islands\" \
                              --trackType CanvasFeatures --trackLabel \"${params.prefix} genomic islands\" \
                              --config '{ "style": { "color": "cyan" }, "displayMode": "compact" }' --out \"data\" ;
  # Remove track for configuration
  [ ! -s all_GIs.bed ] || remove-track.pl --trackLabel \"${params.prefix} genomic islands\" --dir data &> /tmp/error
  # Re-create
  [ ! -s all_GIs.bed ] || echo \' { \"compress\" : 0, \
                                    \"displayMode\" : \"compact\", \
                                    \"key\" : \"Genomic Islands\", \
                                    \"label\" : \"${params.prefix} genomic islands\", \
                                    \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                    \"style\" : { \"className\" : \"feature\", \"color\": \"cyan\" }, \
                                    \"trackType\" : \"CanvasFeatures\", \
                                    \"type\" : \"CanvasFeatures\", \
                                    \"urlTemplate\" : \"tracks/${params.prefix} genomic islands/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json
  # Add track with rRNA sequences
  [ ! -s rrna.gff ] || flatfile-to-json.pl --gff rrna.gff --key \"rRNA Sequences\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} rRNA sequences\" \
  --config '{ "style": { "color": "blue" }, "displayMode": "compact" }' --out \"data\" ;
  remove-track.pl --trackLabel \"${params.prefix} rRNA sequences\" --dir data &> /tmp/error
  [ ! -s rrna.gff ] || echo \' { \"compress\" : 0, \
                                 \"displayMode\" : \"compact\", \
                                 \"key\" : \"rRNA Sequences\", \
                                 \"label\" : \"${params.prefix} rRNA sequences\", \
                                 \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                 \"style\" : { \"className\" : \"feature\", \"color\": \"blue\" }, \
                                 \"trackType\" : \"CanvasFeatures\", \
                                 \"type\" : \"CanvasFeatures\", \
                                 \"urlTemplate\" : \"tracks/${params.prefix} rRNA sequences/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json
  # Add track with efflux features
  [ ! -s efflux ] || flatfile-to-json.pl --gff efflux --key \"Efflux [pumps] features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} efflux features\" --out \"data\" ;
  # Add track with conjugation features
  [ ! -s conjugation ] || flatfile-to-json.pl --gff conjugation --key \"Conjugation related features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} Conjugation related features\" --out \"data\" ;
  # Format bedGraphs
  ## cpg
  [ ! -s cpg ] || bedGraphToBigWig cpg chr.sizes data/cpg.bw ;
  ## gpc
  [ ! -s gpc ] || bedGraphToBigWig gpc chr.sizes data/gpc.bw ;
  ## dam
  [ ! -s dam ] || bedGraphToBigWig dam chr.sizes data/dam.bw ;
  ## dcm
  [ ! -s dcm ] || bedGraphToBigWig dcm chr.sizes data/dcm.bw ;
  # Add BigWigs
  ## cpg
  [ ! -s data/cpg.bw ] || add-bw-track.pl --bw_url cpg.bw --plot --label "CpG Methylations" --key "CpG Methylations" --category "Methylations" \
  --pos_color blue ;
  ## gpc
  [ ! -s data/gpc.bw ] || add-bw-track.pl --bw_url gpc.bw --plot --label "GpC Methylations" --key "GpC Methylations" --category "Methylations" \
  --pos_color purple ;
  ## dam
  [ ! -s data/dam.bw ] || add-bw-track.pl --bw_url dam.bw --plot --label "Dam Methylations" --key "Dam Methylations" --category "Methylations" \
  --pos_color pink ;
  ## dcm
  [ ! -s data/dcm.bw ] || add-bw-track.pl --bw_url dcm.bw --plot --label "Dcm Methylations" --key "Dcm Methylations" --category "Methylations" \
  --pos_color cyan ;
  """
}
