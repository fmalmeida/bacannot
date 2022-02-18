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
include { ISLANDPATH } from '../modules/MGEs/islandpath.nf'
include { DRAW_GIS   } from '../modules/MGEs/draw_gis.nf'

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
include { CUSTOM_DATABASE        } from '../modules/generic/custom_database.nf'
include { CUSTOM_DATABASE_REPORT } from '../modules/generic/custom_database_report.nf'
include { GET_NCBI_PROTEIN       } from '../modules/generic/ncbi_protein.nf'

// Merging annotation in GFF
include { MERGE_ANNOTATIONS } from '../modules/generic/merge_annotations.nf'

// Convert GFF to GBK
include { GFF2GBK } from '../modules/generic/gff2gbk.nf'

// Convert GFF to SQL
include { CREATE_SQL } from '../modules/generic/gff2sql.nf'

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
    dbs_ch
    custom_db
    ncbi_accs
  
  main:

      // generate channel branches
      // now we create the filtered channels
      input_ch.branch{
        unicycler_ch:  it[1] == 'unicycler'
        flye_ch:       it[1] == 'flye'
        annotation_ch: it[1] == "annotation"
      }.set { parsed_inputs }

      // Step 0 -- Run unicycler when necessary
      UNICYCLER( parsed_inputs.unicycler_ch )

      // Step 0 --  Run flye when necessary
      FLYE( parsed_inputs.flye_ch )

      // First step -- Prokka annotation
      PROKKA( 
        parsed_inputs.annotation_ch.mix(FLYE.out[1], UNICYCLER.out[1]), dbs_ch 
      )

      // Second step -- MLST analysis
      MLST( PROKKA.out[3], dbs_ch )

      // Third step -- rRNA annotation
      BARRNAP( PROKKA.out[3] )

      // Fouth step -- calculate GC content for JBrowse
      COMPUTE_GC( PROKKA.out[3] )

      // Fifth step -- run kofamscan
      if (params.skip_kofamscan == false) {
        KOFAMSCAN( PROKKA.out[4], dbs_ch )
        KEGG_DECODER( KOFAMSCAN.out[1] )
        kofamscan_output_ch = KOFAMSCAN.out[1]
        kegg_decoder_svg_ch = KEGG_DECODER.out[1]
      } else {
        kofamscan_output_ch = Channel.empty()
        kegg_decoder_svg_ch = Channel.empty()
      }

      /*
          Sixth step -- MGE, Virulence and AMR annotations
      */

      // plasmids
      if (params.skip_plasmid_search == false) {  
        // plasmidfinder
        PLASMIDFINDER( PROKKA.out[3], dbs_ch )
        plasmidfinder_output_ch = PLASMIDFINDER.out[1]
        // platon
        PLATON( PROKKA.out[3], dbs_ch )
        platon_output_ch = PLATON.out[1]
      } else {
        plasmidfinder_output_ch = Channel.empty()
        platon_output_ch = Channel.empty()
      }

      // IslandPath software
      ISLANDPATH(PROKKA.out[2])

      // Virulence search
      if (params.skip_virulence_search == false) {     
        // VFDB
        VFDB( PROKKA.out[5], dbs_ch )
        vfdb_output_ch = VFDB.out[1]
        // Victors db
        VICTORS( PROKKA.out[4], dbs_ch )
        victors_output_ch = VICTORS.out[1]
      } else {
        vfdb_output_ch    = Channel.empty()
        victors_output_ch = Channel.empty()
      }

      // Prophage search
      if (params.skip_prophage_search == false) {
        // PHAST db
        PHAST( PROKKA.out[4], dbs_ch )
        phast_output_ch = PHAST.out[1]
        // Phigaro software
        PHIGARO( PROKKA.out[3], dbs_ch )
        phigaro_output_1_ch = PHIGARO.out[0]
        phigaro_output_2_ch = PHIGARO.out[1]
        // PhiSpy
        PHISPY( PROKKA.out[2] )
        phispy_output_ch = PHISPY.out[1]
      } else {
        phast_output_ch     = Channel.empty()
        phigaro_output_1_ch = Channel.empty()
        phigaro_output_2_ch = Channel.empty()
        phispy_output_ch    = Channel.empty()
      }

      // ICEs search
      if (params.skip_iceberg_search == false) {
        // ICEberg db
        ICEBERG( PROKKA.out[4], PROKKA.out[3], dbs_ch )
        iceberg_output_ch   = ICEBERG.out[1]
        iceberg_output_2_ch = ICEBERG.out[2]
      } else {
        iceberg_output_ch   = Channel.empty()
        iceberg_output_2_ch = Channel.empty()
      }

      // AMR search
      if (params.skip_resistance_search == false) {
        // AMRFinderPlus
        AMRFINDER( PROKKA.out[4], dbs_ch )
        amrfinder_output_ch = AMRFINDER.out[0]
        // CARD-RGI
        CARD_RGI( PROKKA.out[4], dbs_ch )
        rgi_output_ch        = CARD_RGI.out[2]
        rgi_output_parsed_ch = CARD_RGI.out[1]
        rgi_heatmap_ch       = CARD_RGI.out[3]
        // ARGMiner
        ARGMINER( PROKKA.out[4], dbs_ch )
        argminer_output_ch = ARGMINER.out[0]
        // Resfinder
        RESFINDER( PROKKA.out[7], dbs_ch )
        resfinder_output_1_ch   = RESFINDER.out[0]
        resfinder_output_2_ch   = RESFINDER.out[1]
        resfinder_phenotable_ch = RESFINDER.out[2]
        resfinder_gff_ch        = RESFINDER.out[3]
      } else {
        rgi_output_ch           = Channel.empty()
        rgi_output_parsed_ch    = Channel.empty()
        rgi_heatmap_ch          = Channel.empty()
        amrfinder_output_ch     = Channel.empty()
        argminer_output_ch      = Channel.empty()
        resfinder_output_1_ch   = Channel.empty()
        resfinder_output_2_ch   = Channel.empty()
        resfinder_phenotable_ch = Channel.empty()
        resfinder_gff_ch        = Channel.empty()
      }

      /*
          Seventh step -- Methylation call
      */
      CALL_METHYLATION( PROKKA.out[6] )
      methylation_out_1_ch = CALL_METHYLATION.out[2]
      methylation_out_2_ch = CALL_METHYLATION.out[3]

      /*

          Additional steps created after main releases

       */

      // species identification
      REFSEQ_MASHER( PROKKA.out[3] )

      // IS identification
      DIGIS( PROKKA.out[3].join(PROKKA.out[2]) )

      // antiSMASH
      if (params.skip_antismash == false) {
        ANTISMASH( PROKKA.out[2], dbs_ch )
        antismash_output_ch = ANTISMASH.out[0]
      } else {
        antismash_output_ch = Channel.empty()
      }

      // sequenceserver
      SEQUENCESERVER(
        PROKKA.out[3].join(PROKKA.out[5])
                     .join(PROKKA.out[4])
      )

      // custom databases annotation
      ch_custom_databases_annotations = Channel.empty()
      if (params.custom_db || params.ncbi_proteins) {
        GET_NCBI_PROTEIN( ncbi_accs )
        CUSTOM_DATABASE(
          PROKKA.out[1].join(PROKKA.out[3]),
          custom_db.mix(GET_NCBI_PROTEIN.out[0])
        )
        ch_custom_databases_annotations = CUSTOM_DATABASE.out[1].groupTuple()
      }

      /*
          Eighth step -- Merge all annotations
      */
      // prefix is 0
      // annotations_files_ch = PROKKA.out[3] // 1
      //                        .join(PROKKA.out[1]) // 2
      //                        .join(MLST.out[0]) // 3
      //                        .join(BARRNAP.out[0]) // 4
      //                        .join(COMPUTE_GC.out[0]) // 5
      //                        .join(kofamscan_output_ch,             remainder: true) // 6
      //                        .join(vfdb_output_ch,                  remainder: true) // 7
      //                        .join(victors_output_ch,               remainder: true) // 8
      //                        .join(amrfinder_output_ch,             remainder: true) // 9
      //                        .join(resfinder_gff_ch,                remainder: true) // 10
      //                        .join(rgi_output_ch,                   remainder: true) // 11
      //                        .join(iceberg_output_ch,               remainder: true) // 12
      //                        .join(phast_output_ch,                 remainder: true) // 13
      //                        .join(phigaro_output_2_ch,             remainder: true) // 14
      //                        .join(ISLANDPATH.out[0],               remainder: true) // 15
      //                        .join(DIGIS.out[1],                    remainder: true) // 16
      //                        .join(ISLANDPATH.out[0],               remainder: true) // 17
      //                        .join(ch_custom_databases_annotations, remainder: true) // 18
      //                        .toList() // transforms in single list with everything

      // Contatenation of annotations in a single GFF file
      MERGE_ANNOTATIONS( 
        PROKKA.out[1].join(kofamscan_output_ch,             remainder: true)
                     .join(vfdb_output_ch,                  remainder: true)
                     .join(victors_output_ch,               remainder: true)
                     .join(amrfinder_output_ch,             remainder: true)
                     .join(resfinder_gff_ch,                remainder: true)
                     .join(rgi_output_ch,                   remainder: true)
                     .join(iceberg_output_ch,               remainder: true)
                     .join(phast_output_ch,                 remainder: true)
                     .join(DIGIS.out[1],                    remainder: true)
                     .join(ch_custom_databases_annotations, remainder: true) 
      )

      // Plot genomic islands
      DRAW_GIS( MERGE_ANNOTATIONS.out[0] )

      // Convert GFF file to GBK file
      GFF2GBK( MERGE_ANNOTATIONS.out[0].join(PROKKA.out[3]) )

      // Convert GFF file to sqldb
      CREATE_SQL(
        MERGE_ANNOTATIONS.out[0].join(PROKKA.out[5])
                                .join(PROKKA.out[4])
                                .join(PROKKA.out[3])
                                .join(DIGIS.out[2] )
      )

      // /*
      //     Final step -- Create genome browser and reports
      // */

      // // Grab inputs needed for JBrowse step
      // jbrowse_input_ch = MERGE_ANNOTATIONS.out[0].join(annotations_files_ch, remainder: true)
      //                                         .join(methylation_out_1_ch, remainder: true)
      //                                         .join(methylation_out_2_ch, remainder: true)
      //                                         .join(phispy_output_ch,     remainder: true)
      //                                         .join(MERGE_ANNOTATIONS.out[8], remainder: true) // parsed and changed digIS
      //                                         .join(antismash_output_ch,  remainder: true)
      // // Jbrowse Creation
      // JBROWSE(jbrowse_input_ch)

      // // Render reports
      // CUSTOM_DATABASE_REPORT( CUSTOM_DATABASE.out[0].join(CUSTOM_DATABASE.out[1]) )
      // REPORT(jbrowse_input_ch.join(rgi_output_parsed_ch,    remainder: true)
      //                     .join(rgi_heatmap_ch,          remainder: true)
      //                     .join(argminer_output_ch,      remainder: true)
      //                     .join(iceberg_output_2_ch,     remainder: true)
      //                     .join(plasmidfinder_output_ch, remainder: true)
      //                     .join(resfinder_output_1_ch,   remainder: true)
      //                     .join(resfinder_output_2_ch,   remainder: true)
      //                     .join(resfinder_phenotable_ch, remainder: true)
      //                     .join(DRAW_GIS.out[1],      remainder: true)
      //                     .join(phigaro_output_1_ch,     remainder: true)
      //                     .join(platon_output_ch,        remainder: true)
      //                     .join(PROKKA.out[8],        remainder: true)
      //                     .join(kegg_decoder_svg_ch,     remainder: true)
      //                     .join(REFSEQ_MASHER.out[0], remainder: true))

}
