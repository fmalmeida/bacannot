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

// Bakta annotation
include { BAKTA } from '../modules/generic/bakta.nf'

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

      // First step -- Prokka or Bakta annotation
      if (params.bakta_db) {
        BAKTA(
          parsed_inputs.annotation_ch.mix(FLYE.out[1], UNICYCLER.out[1]),
          file(params.bakta_db, checkIfExists: true)  
        )
        annotation_out_ch = BAKTA.out
      } else {
        PROKKA(parsed_inputs.annotation_ch.mix(FLYE.out[1], UNICYCLER.out[1]), dbs_ch)
        annotation_out_ch = PROKKA.out
      }

      // Second step -- MLST analysis
      MLST( annotation_out_ch[3], dbs_ch )

      // Third step -- rRNA annotation
      BARRNAP( annotation_out_ch[3] )

      // Fouth step -- calculate GC content for JBrowse
      COMPUTE_GC( annotation_out_ch[3] )

      // Fifth step -- run kofamscan
      if (params.skip_kofamscan == false) {
        KOFAMSCAN( annotation_out_ch[4], dbs_ch )
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
        PLASMIDFINDER( annotation_out_ch[3], dbs_ch )
        plasmidfinder_output_ch = PLASMIDFINDER.out[1]
        // platon
        PLATON( annotation_out_ch[3], dbs_ch )
        platon_output_ch = PLATON.out[1]
      } else {
        plasmidfinder_output_ch = Channel.empty()
        platon_output_ch = Channel.empty()
      }

      // IslandPath software
      ISLANDPATH(annotation_out_ch[2])

      // Virulence search
      if (params.skip_virulence_search == false) {     
        // VFDB
        VFDB( annotation_out_ch[5], dbs_ch )
        vfdb_output_ch = VFDB.out[1]
        // Victors db
        VICTORS( annotation_out_ch[4], dbs_ch )
        victors_output_ch = VICTORS.out[1]
      } else {
        vfdb_output_ch    = Channel.empty()
        victors_output_ch = Channel.empty()
      }

      // Prophage search
      if (params.skip_prophage_search == false) {
        // PHAST db
        PHAST( annotation_out_ch[4], dbs_ch )
        phast_output_ch = PHAST.out[1]
        // Phigaro software
        PHIGARO( annotation_out_ch[3], dbs_ch )
        phigaro_output_tsv_ch = PHIGARO.out[0]
        phigaro_output_bed_ch = PHIGARO.out[1]
        // PhiSpy
        PHISPY( annotation_out_ch[2] )
        phispy_output_ch = PHISPY.out[1]
      } else {
        phast_output_ch       = Channel.empty()
        phigaro_output_tsv_ch = Channel.empty()
        phigaro_output_bed_ch = Channel.empty()
        phispy_output_ch      = Channel.empty()
      }

      // ICEs search
      if (params.skip_iceberg_search == false) {
        // ICEberg db
        ICEBERG( annotation_out_ch[4], annotation_out_ch[3], dbs_ch )
        iceberg_output_blastp_ch   = ICEBERG.out[1]
        iceberg_output_blastn_ch   = ICEBERG.out[2]
      } else {
        iceberg_output_blastp_ch   = Channel.empty()
        iceberg_output_blastn_ch   = Channel.empty()
      }

      // AMR search
      if (params.skip_resistance_search == false) {
        // AMRFinderPlus
        AMRFINDER( annotation_out_ch[4], dbs_ch )
        amrfinder_output_ch = AMRFINDER.out[0]
        // CARD-RGI
        CARD_RGI( annotation_out_ch[4], dbs_ch )
        rgi_output_ch        = CARD_RGI.out[2]
        rgi_output_parsed_ch = CARD_RGI.out[1]
        rgi_heatmap_ch       = CARD_RGI.out[3]
        // ARGMiner
        ARGMINER( annotation_out_ch[4], dbs_ch )
        argminer_output_ch = ARGMINER.out[0]
        // Resfinder
        RESFINDER( annotation_out_ch[7], dbs_ch )
        resfinder_output_tab_ch           = RESFINDER.out[0]
        resfinder_output_pointfinder_ch   = RESFINDER.out[1]
        resfinder_phenotable_ch           = RESFINDER.out[2]
        resfinder_gff_ch                  = RESFINDER.out[3]
      } else {
        rgi_output_ch                     = Channel.empty()
        rgi_output_parsed_ch              = Channel.empty()
        rgi_heatmap_ch                    = Channel.empty()
        amrfinder_output_ch               = Channel.empty()
        argminer_output_ch                = Channel.empty()
        resfinder_output_tab_ch           = Channel.empty()
        resfinder_output_pointfinder_ch   = Channel.empty()
        resfinder_phenotable_ch           = Channel.empty()
        resfinder_gff_ch                  = Channel.empty()
      }

      /*
          Seventh step -- Methylation call
      */
      CALL_METHYLATION( annotation_out_ch[6] )
      methylation_out_1_ch = CALL_METHYLATION.out[2]
      methylation_out_2_ch = CALL_METHYLATION.out[3]

      /*

          Additional steps created after main releases

       */

      // species identification
      REFSEQ_MASHER( annotation_out_ch[3] )

      // IS identification
      DIGIS( annotation_out_ch[3].join(annotation_out_ch[2]) )

      // antiSMASH
      if (params.skip_antismash == false) {
        ANTISMASH( annotation_out_ch[2], dbs_ch )
        antismash_output_ch = ANTISMASH.out[0]
      } else {
        antismash_output_ch = Channel.empty()
      }

      // sequenceserver
      SEQUENCESERVER(
        annotation_out_ch[3].join(annotation_out_ch[5])
                     .join(annotation_out_ch[4])
      )

      // custom databases annotation
      ch_custom_databases_annotations = Channel.empty()
      if (params.custom_db || params.ncbi_proteins) {
        GET_NCBI_PROTEIN( ncbi_accs )
        CUSTOM_DATABASE(
          annotation_out_ch[1].join(annotation_out_ch[3]),
          custom_db.mix(GET_NCBI_PROTEIN.out[0])
        )
        ch_custom_databases_annotations = CUSTOM_DATABASE.out[1].groupTuple()
      }

      /*
          Eighth step -- Merge all annotations
      */
      MERGE_ANNOTATIONS( 
        annotation_out_ch[1].join(kofamscan_output_ch,             remainder: true)
                     .join(vfdb_output_ch,                  remainder: true)
                     .join(victors_output_ch,               remainder: true)
                     .join(amrfinder_output_ch,             remainder: true)
                     .join(resfinder_gff_ch,                remainder: true)
                     .join(rgi_output_ch,                   remainder: true)
                     .join(iceberg_output_blastp_ch,        remainder: true)
                     .join(phast_output_ch,                 remainder: true)
                     .join(DIGIS.out[1],                    remainder: true)
                     .join(ch_custom_databases_annotations, remainder: true) 
      )

      /*
          Final step -- Create genome browser and reports' files
      */
      // Plot genomic islands
      DRAW_GIS( MERGE_ANNOTATIONS.out[0].join(ISLANDPATH.out[0]) )

      // Convert GFF file to GBK file
      GFF2GBK( MERGE_ANNOTATIONS.out[0].join(annotation_out_ch[3]) )

      // Convert GFF file to sqldb
      CREATE_SQL(
        MERGE_ANNOTATIONS.out[0].join(annotation_out_ch[5])
                                .join(annotation_out_ch[4])
                                .join(annotation_out_ch[3])
                                .join(DIGIS.out[2] )
      )

      JBROWSE(
        MERGE_ANNOTATIONS.out[0].join(annotation_out_ch[3])
                                .join(annotation_out_ch[1])
                                .join(BARRNAP.out[0])
                                .join(COMPUTE_GC.out[0])
                                .join(resfinder_gff_ch,     remainder: true)
                                .join(phigaro_output_bed_ch,remainder: true)
                                .join(ISLANDPATH.out[0],    remainder: true)
                                .join(methylation_out_1_ch, remainder: true)
                                .join(methylation_out_2_ch, remainder: true)
                                .join(phispy_output_ch,     remainder: true)
                                .join(MERGE_ANNOTATIONS.out[1]) // parsed digIS
                                .join(antismash_output_ch,  remainder: true)
                                .join(MERGE_ANNOTATIONS.out[2].groupTuple(), remainder: true) // parsed custom db
      )

      // Render reports
      if (params.custom_db || params.ncbi_proteins) {
        CUSTOM_DATABASE_REPORT( CUSTOM_DATABASE.out[0].join( MERGE_ANNOTATIONS.out[0], remainder:true ) )
      }
      REPORT(
        annotation_out_ch[8].join(MERGE_ANNOTATIONS.out[0])
                     .join(BARRNAP.out[0])
                     .join(MLST.out[0])
                     .join(kegg_decoder_svg_ch,             remainder: true)
                     .join(REFSEQ_MASHER.out[0])
                     .join(amrfinder_output_ch,             remainder: true)
                     .join(rgi_output_ch,                   remainder: true)
                     .join(rgi_output_parsed_ch,            remainder: true)
                     .join(rgi_heatmap_ch,                  remainder: true)
                     .join(argminer_output_ch,              remainder: true)
                     .join(resfinder_output_tab_ch,         remainder: true)
                     .join(resfinder_output_pointfinder_ch, remainder: true)
                     .join(resfinder_phenotable_ch,         remainder: true)
                     .join(vfdb_output_ch,                  remainder: true)
                     .join(victors_output_ch,               remainder: true)
                     .join(phigaro_output_tsv_ch,           remainder: true)
                     .join(phispy_output_ch,                remainder: true)
                     .join(iceberg_output_blastp_ch,        remainder: true)
                     .join(iceberg_output_blastn_ch,        remainder: true)
                     .join(plasmidfinder_output_ch,         remainder: true)
                     .join(platon_output_ch,                remainder: true)
                     .join(DRAW_GIS.out[1],                 remainder: true)
                     .join(phast_output_ch,                 remainder: true)
                     .join(MERGE_ANNOTATIONS.out[1]) // parsed digIS
      )

}
