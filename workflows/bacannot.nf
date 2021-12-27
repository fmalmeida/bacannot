/*
 * Include modules (Execution setup)
 */

// Unicycler assembly
include { UNICYCLER } from '../modules/assembly/unicycler.nf'

// Flye assembly
include { FLYE } from '../modules/assembly/flye.nf'

// Species identification
include { REFSEQ_MASHER } from '../modules/generic/mash.nf'

// Prokka annotation
include { PROKKA } from '../modules/generic/prokka.nf'

// MLST annotation
include { MLST } from '../modules/generic/mlst.nf'

// rRNA annotation
include { BARRNAP } from '../modules/generic/barrnap.nf'

// Calculate GC content
include { COMPUTE_GC } from '../modules/generic/compute_gc.nf'

// KOFAM annotation
include { KOFAMSCAN } from '../modules/KOs/kofamscan.nf'

// KEGG decoder
include { KEGG_DECODER } from '../modules/KOs/kegg-decoder.nf'

// Plasmid annotation with plasmidfinder
include { PLASMIDFINDER } from '../modules/MGEs/plasmidfinder.nf'

// Plasmid annotation with platon
include { PLATON } from '../modules/MGEs/platon.nf'

// Virulence annotation with VFDB
include { VFDB } from '../modules/virulence/vfdb.nf'

// Virulence annotation with Victors
include { VICTORS } from '../modules/virulence/victors.nf'

// Prophage annotation with PHAST
include { PHAST } from '../modules/prophages/phast.nf'

// Prophage annotation with PHIGARO
include { PHIGARO } from '../modules/prophages/phigaro.nf'

// Prophage annotation with phispy
include { PHISPY } from '../modules/prophages/phispy.nf'

// ICE annotation with ICEberg db
include { ICEBERG } from '../modules/MGEs/iceberg.nf'

// Genomic Islands detection with Islandpath-DIMOB
include { FIND_GIS } from '../modules/MGEs/islandPath_DIMOB.nf'
include { DRAW_GIS } from '../modules/MGEs/draw_gis.nf'

// IS identification
include { DIGIS } from '../modules/MGEs/digIS.nf'

// AMR annotation with ARGMiner
include { ARGMINER } from '../modules/resistance/argminer.nf'

// AMR annotation with Resfinder
include { RESFINDER } from '../modules/resistance/resfinder.nf'

// AMR annotation with AMRFinderPlus
include { AMRFINDER } from '../modules/resistance/amrfinder.nf'

// AMR annotation with CARD-RGI
include { CARD_RGI } from '../modules/resistance/rgi_annotation.nf'

// Methylation calling (Nanopolish)
include { CALL_METHYLATION } from '../modules/generic/methylation.nf'

// User's custom db annotation
include { CUSTOM_BLAST } from '../modules/generic/custom_blast.nf'
include { CUSTOM_BLAST_REPORT } from '../modules/generic/custom_blast_report.nf'

// Merging annotation in GFF
include { MERGE_ANNOTATIONS } from '../modules/generic/merge_annotations.nf'

// Convert GFF to GBK
include { GFF2GBK } from '../modules/generic/gff2gbk.nf'

// Convert GFF to SQL
include { CREATE_SQL } from '../modules/generic/gff2sql.nf'

// Bedtools gff merge
include { GFF_MERGE } from '../modules/generic/merge_gff.nf'

// JBrowse
include { JBROWSE } from '../modules/generic/jbrowse.nf'

// Output reports
include { REPORT } from '../modules/generic/reports.nf'

// sequenceserver generation
include { SEQUENCESERVER } from '../modules/generic/sequenceserver.nf'

// antiSMASH
include { ANTISMASH } from '../modules/generic/antismash.nf'

/*
    DEF WORKFLOW
*/

workflow BACANNOT {
  take:
    input_ch
    custom_db
  
  main:

      // generate channel branches
      // now we create the filtered channels
      input_ch.branch{
        unicycler_ch:  it[1] == 'unicycler'
        flye_ch:       it[1] == 'flye'
        annotation_ch: it[1] == "annotation"
      }.set { parsed_inputs }

      // Step 0 -- Run unicycler when necessary
      UNICYCLER(parsed_inputs.unicycler_ch)

      // Step 0 --  Run flye when necessary
      FLYE(parsed_inputs.flye_ch)

      // First step -- Prokka annotation
      PROKKA(parsed_inputs.annotation_ch.mix(FLYE.out[1], UNICYCLER.out[1]))

      // Second step -- MLST analysis
      MLST(PROKKA.out[3])

      // Third step -- rRNA annotation
      BARRNAP(PROKKA.out[3])

      // Fouth step -- calculate GC content for JBrowse
      COMPUTE_GC(PROKKA.out[3])

      // Fifth step -- run kofamscan
      if (params.skip_kofamscan == false) {
        KOFAMSCAN(PROKKA.out[4])
        KEGG_DECODER(KOFAMSCAN.out[1])
        kofamscan_output_ch = KOFAMSCAN.out[1]
        kegg_decoder_svg_ch = KEGG_DECODER.out[1]
      } else {
        kofamscan_output_ch = Channel.empty()
        kegg_decoder_svg_ch = Channel.empty()
      }

      /*
          Sixth step -- MGE, Virulence and AMR annotations
      */

      // Plasmid finder
      if (params.skip_plasmid_search == false) {
        // plasmidfinder
        PLASMIDFINDER(PROKKA.out[3])
        plasmidfinder_output_ch = PLASMIDFINDER.out[1]
        // platon
        PLATON(PROKKA.out[3])
        platon_output_ch = PLATON.out[1]
      } else {
        plasmidfinder_output_ch = Channel.empty()
        platon_output_ch = Channel.empty()
      }

      // IslandPath software
      FIND_GIS(PROKKA.out[2])

      // Virulence search
      if (params.skip_virulence_search == false) {
        // VFDB
        VFDB(PROKKA.out[5])
        vfdb_output_ch = VFDB.out[1]
        // Victors db
        VICTORS(PROKKA.out[4])
        victors_output_ch = VICTORS.out[1]
      } else {
        vfdb_output_ch = Channel.empty()
        victors_output_ch = Channel.empty()
      }

      // Prophage search
      if (params.skip_prophage_search == false) {
        // PHAST db
        PHAST(PROKKA.out[4])
        phast_output_ch = PHAST.out[1]
        // Phigaro software
        PHIGARO(PROKKA.out[3])
        phigaro_output_1_ch = PHIGARO.out[0]
        phigaro_output_2_ch = PHIGARO.out[1]
        // PhiSpy
        PHISPY(PROKKA.out[2])
        phispy_output_ch = PHISPY.out[1]
      } else {
        phast_output_ch = Channel.empty()
        phigaro_output_1_ch = Channel.empty()
        phigaro_output_2_ch = Channel.empty()
        phispy_output_ch = Channel.empty()
      }

      // ICEs search
      if (params.skip_iceberg_search == false) {
        // ICEberg db
        ICEBERG(PROKKA.out[4], PROKKA.out[3])
        iceberg_output_ch = ICEBERG.out[1]
        iceberg_output_2_ch = ICEBERG.out[2]
      } else {
        iceberg_output_ch = Channel.empty()
        iceberg_output_2_ch = Channel.empty()
      }

      // AMR search
      if (params.skip_resistance_search == false) {
        // AMRFinderPlus
        AMRFINDER(PROKKA.out[4])
        amrfinder_output_ch = AMRFINDER.out[0]
        // CARD-RGI
        CARD_RGI(PROKKA.out[4])
        rgi_output_ch = CARD_RGI.out[2]
        rgi_output_parsed_ch = CARD_RGI.out[1]
        rgi_heatmap_ch = CARD_RGI.out[3]
        // ARGMiner
        ARGMINER(PROKKA.out[4])
        argminer_output_ch = ARGMINER.out[0]
        // Resfinder
        RESFINDER(PROKKA.out[7])
        resfinder_output_1_ch = RESFINDER.out[0]
        resfinder_output_2_ch = RESFINDER.out[1]
        resfinder_phenotable_ch = RESFINDER.out[2]
        resfinder_gff_ch = RESFINDER.out[3]
      } else {
        rgi_output_ch = Channel.empty()
        rgi_output_parsed_ch = Channel.empty()
        rgi_heatmap_ch = Channel.empty()
        amrfinder_output_ch = Channel.empty()
        argminer_output_ch = Channel.empty()
        resfinder_output_1_ch = Channel.empty()
        resfinder_output_2_ch = Channel.empty()
        resfinder_phenotable_ch = Channel.empty()
        resfinder_gff_ch = Channel.empty()
      }

      /*
          Seventh step -- Methylation call
      */
      CALL_METHYLATION(PROKKA.out[6])
      methylation_out_1_ch = CALL_METHYLATION.out[2]
      methylation_out_2_ch = CALL_METHYLATION.out[3]

      /*

          Additional steps created after main releases

       */

      // species identification
      REFSEQ_MASHER(PROKKA.out[3])

      // IS identification
      DIGIS(PROKKA.out[3].join(PROKKA.out[2]))

      // antiSMASH
      if (params.skip_antismash == false) {
        ANTISMASH(PROKKA.out[2])
        antismash_output_ch = ANTISMASH.out[0]
      } else {
        antismash_output_ch = Channel.empty()
      }

      // sequenceserver
      SEQUENCESERVER(PROKKA.out[3].join(PROKKA.out[5]).join(PROKKA.out[4]))

      /*
          Eighth step -- Merge all annotations with the same Prefix value in a single Channel
      */
      annotations_files_ch = PROKKA.out[3].join(PROKKA.out[1])
                                       .join(MLST.out[0])
                                       .join(BARRNAP.out[0])
                                       .join(COMPUTE_GC.out[0])
                                       .join(kofamscan_output_ch, remainder: true)
                                       .join(vfdb_output_ch,      remainder: true)
                                       .join(victors_output_ch,   remainder: true)
                                       .join(amrfinder_output_ch, remainder: true)
                                       .join(resfinder_gff_ch,    remainder: true)
                                       .join(rgi_output_ch,       remainder: true)
                                       .join(iceberg_output_ch,   remainder: true)
                                       .join(phast_output_ch,     remainder: true)
                                       .join(phigaro_output_2_ch, remainder: true)
                                       .join(FIND_GIS.out[0],  remainder: true)

      // Contatenation of annotations in a single GFF file
      MERGE_ANNOTATIONS(ANNOTATIONS_FILES.join(DIGIS.out[1],     remainder: true))

      // Plot genomic islands
      DRAW_GIS(MERGE_ANNOTATIONS.out[0].join(FIND_GIS.out[0]))

      // Convert GFF file to GBK file
      GFF2GBK(MERGE_ANNOTATIONS.out[0].join(PROKKA.out[3]))

      // Convert GFF file to sqldb
      CREATE_SQL(MERGE_ANNOTATIONS.out[0].join(PROKKA.out[5])
                                         .join(PROKKA.out[4])
                                         .join(PROKKA.out[3])
                                         .join(DIGIS.out[2]))

      // User wants to merge the final gff file?
      if (params.bedtools_merge_distance) {
        GFF_MERGE(MERGE_ANNOTATIONS.out[0])
      }

      /*

          Nineth step -- Perform users custom annotation

      */
      if (params.custom_db) {
        CUSTOM_BLAST(MERGE_ANNOTATIONS.out[0].join(PROKKA.out[3]), custom_db)
        CUSTOM_BLAST_REPORT(CUSTOM_BLAST.out[0])
      }

      /*
          Final step -- Create genome browser and reports
      */

      // Grab inputs needed for JBrowse step
      jbrowse_input_ch = merge_annotations.out[0].join(annotations_files_ch, remainder: true)
                                              .join(methylation_out_1_ch, remainder: true)
                                              .join(methylation_out_2_ch, remainder: true)
                                              .join(phispy_output_ch,     remainder: true)
                                              .join(MERGE_ANNOTATIONS.out[8], remainder: true) // parsed and changed digIS
                                              .join(antismash_output_ch,  remainder: true)
      // Jbrowse Creation
      JBROWSE(jbrowse_input_ch)

      // Render reports
      REPORT(jbrowse_input_ch.join(rgi_output_parsed_ch,    remainder: true)
                          .join(rgi_heatmap_ch,          remainder: true)
                          .join(argminer_output_ch,      remainder: true)
                          .join(iceberg_output_2_ch,     remainder: true)
                          .join(plasmidfinder_output_ch, remainder: true)
                          .join(resfinder_output_1_ch,   remainder: true)
                          .join(resfinder_output_2_ch,   remainder: true)
                          .join(resfinder_phenotable_ch, remainder: true)
                          .join(DRAW_GIS.out[1],      remainder: true)
                          .join(phigaro_output_1_ch,     remainder: true)
                          .join(platon_output_ch,        remainder: true)
                          .join(PROKKA.out[8],        remainder: true)
                          .join(kegg_decoder_svg_ch,     remainder: true)
                          .join(REFSEQ_MASHER.out[0], remainder: true))

}
