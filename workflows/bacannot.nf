/*
 * Include modules (Execution setup)
 */

include { UNICYCLER              } from '../modules/assembly/unicycler.nf'
include { FLYE                   } from '../modules/assembly/flye.nf'
include { REFSEQ_MASHER          } from '../modules/generic/mash.nf'
include { PROKKA                 } from '../modules/generic/prokka.nf'
include { BAKTA                  } from '../modules/generic/bakta.nf'
include { MLST                   } from '../modules/generic/mlst.nf'
include { BARRNAP                } from '../modules/generic/barrnap.nf'
include { COMPUTE_GC             } from '../modules/generic/compute_gc.nf'
include { KOFAMSCAN              } from '../modules/KOs/kofamscan.nf'
include { KEGG_DECODER           } from '../modules/KOs/kegg-decoder.nf'
include { PLASMIDFINDER          } from '../modules/MGEs/plasmidfinder.nf'
include { PLATON                 } from '../modules/MGEs/platon.nf'
include { VFDB                   } from '../modules/virulence/vfdb.nf'
include { VICTORS                } from '../modules/virulence/victors.nf'
include { PHAST                  } from '../modules/prophages/phast.nf'
include { PHIGARO                } from '../modules/prophages/phigaro.nf'
include { PHISPY                 } from '../modules/prophages/phispy.nf'
include { ICEBERG                } from '../modules/MGEs/iceberg.nf'
include { ISLANDPATH             } from '../modules/MGEs/islandpath.nf'
include { DRAW_GIS               } from '../modules/MGEs/draw_gis.nf'
include { DIGIS                  } from '../modules/MGEs/digIS.nf'
include { ARGMINER               } from '../modules/resistance/argminer.nf'
include { RESFINDER              } from '../modules/resistance/resfinder.nf'
include { AMRFINDER              } from '../modules/resistance/amrfinder.nf'
include { CARD_RGI               } from '../modules/resistance/rgi_annotation.nf'
include { CALL_METHYLATION       } from '../modules/generic/methylation.nf'
include { CUSTOM_DATABASE        } from '../modules/generic/custom_database.nf'
include { CUSTOM_DATABASE_REPORT } from '../modules/generic/custom_database_report.nf'
include { GET_NCBI_PROTEIN       } from '../modules/generic/ncbi_protein.nf'
include { MERGE_ANNOTATIONS      } from '../modules/generic/merge_annotations.nf'
include { GFF2GBK                } from '../modules/generic/gff2gbk.nf'
include { CREATE_SQL             } from '../modules/generic/gff2sql.nf'
include { JBROWSE                } from '../modules/generic/jbrowse.nf'
include { REPORT                 } from '../modules/generic/reports.nf'
include { SEQUENCESERVER         } from '../modules/generic/sequenceserver.nf'
include { ANTISMASH              } from '../modules/generic/antismash.nf'

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
          parsed_inputs.annotation_ch.mix( FLYE.out.results, UNICYCLER.out.results ),
          file(params.bakta_db, checkIfExists: true)  
        )
        annotation_out_ch = BAKTA.out
      } else {
        PROKKA(
          parsed_inputs.annotation_ch.mix( FLYE.out.results, UNICYCLER.out.results ), 
          dbs_ch
        )
        annotation_out_ch = PROKKA.out
      }

      // Second step -- MLST analysis
      MLST( annotation_out_ch.genome, dbs_ch )

      // Third step -- rRNA annotation
      BARRNAP( annotation_out_ch.genome )

      // Fouth step -- calculate GC content for JBrowse
      COMPUTE_GC( annotation_out_ch.genome )

      // Fifth step -- run kofamscan
      if (params.skip_kofamscan == false) {
        KOFAMSCAN( annotation_out_ch.proteins, dbs_ch )
        KEGG_DECODER( KOFAMSCAN.out.results )
        kofamscan_output_ch = KOFAMSCAN.out.results
        kegg_decoder_svg_ch = KEGG_DECODER.out.results
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
        PLASMIDFINDER( annotation_out_ch.genome, dbs_ch )
        plasmidfinder_output_ch = PLASMIDFINDER.out.results
        // platon
        PLATON( annotation_out_ch.genome, dbs_ch )
        platon_output_ch = PLATON.out.results
      } else {
        plasmidfinder_output_ch = Channel.empty()
        platon_output_ch        = Channel.empty()
      }

      // IslandPath software
      ISLANDPATH( annotation_out_ch.gbk )

      // Virulence search
      if (params.skip_virulence_search == false) {     
        // VFDB
        VFDB( annotation_out_ch.genes, dbs_ch )
        vfdb_output_ch = VFDB.out.results
        // Victors db
        VICTORS( annotation_out_ch.proteins, dbs_ch )
        victors_output_ch = VICTORS.out.results
      } else {
        vfdb_output_ch    = Channel.empty()
        victors_output_ch = Channel.empty()
      }

      // Prophage search
      if (params.skip_prophage_search == false) {
        // PHAST db
        PHAST( annotation_out_ch.proteins, dbs_ch )
        phast_output_ch = PHAST.out.results
        // Phigaro software
        PHIGARO( annotation_out_ch.genome, dbs_ch )
        phigaro_output_tsv_ch = PHIGARO.out.tsv
        phigaro_output_bed_ch = PHIGARO.out.bed
        // PhiSpy
        PHISPY( annotation_out_ch.gbk )
        phispy_output_ch = PHISPY.out.results
      } else {
        phast_output_ch       = Channel.empty()
        phigaro_output_tsv_ch = Channel.empty()
        phigaro_output_bed_ch = Channel.empty()
        phispy_output_ch      = Channel.empty()
      }

      // ICEs search
      if (params.skip_iceberg_search == false) {
        // ICEberg db
        ICEBERG( annotation_out_ch.proteins, annotation_out_ch.genome, dbs_ch )
        iceberg_output_blastp_ch   = ICEBERG.out.results
        iceberg_output_blastn_ch   = ICEBERG.out.genome_summary
      } else {
        iceberg_output_blastp_ch   = Channel.empty()
        iceberg_output_blastn_ch   = Channel.empty()
      }

      // AMR search
      if (params.skip_resistance_search == false) {
        // AMRFinderPlus
        AMRFINDER( annotation_out_ch.proteins, dbs_ch )
        amrfinder_output_ch = AMRFINDER.out.resistance_results
        // CARD-RGI
        CARD_RGI( annotation_out_ch.proteins, dbs_ch )
        rgi_output_ch        = CARD_RGI.out.raw_hits
        rgi_output_parsed_ch = CARD_RGI.out.parsed_hits
        rgi_heatmap_ch       = CARD_RGI.out.heatmap_png
        // ARGMiner
        ARGMINER( annotation_out_ch.proteins, dbs_ch )
        argminer_output_ch = ARGMINER.out.summary
        // Resfinder
        RESFINDER( annotation_out_ch.genome_with_species, dbs_ch )
        resfinder_output_tab_ch           = RESFINDER.out.results
        resfinder_output_pointfinder_ch   = RESFINDER.out.pointfinder_results
        resfinder_phenotable_ch           = RESFINDER.out.pheno_table
        resfinder_gff_ch                  = RESFINDER.out.gff
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
      CALL_METHYLATION( annotation_out_ch.genome_with_fast5 )

      /*

          Additional steps created after main releases

       */

      // species identification
      REFSEQ_MASHER( annotation_out_ch.genome )

      // IS identification
      DIGIS( annotation_out_ch.genome.join(annotation_out_ch.gbk) )

      // antiSMASH
      if (params.skip_antismash == false) {
        ANTISMASH( annotation_out_ch.gbk, dbs_ch )
        antismash_output_ch = ANTISMASH.out.gff
      } else {
        antismash_output_ch = Channel.empty()
      }

      // sequenceserver
      SEQUENCESERVER(
        annotation_out_ch.genome
          .join( annotation_out_ch.genes    )
          .join( annotation_out_ch.proteins )
      )

      // custom databases annotation
      ch_custom_databases_annotations = Channel.empty()
      if (params.custom_db || params.ncbi_proteins) {
        GET_NCBI_PROTEIN( ncbi_accs )
        CUSTOM_DATABASE(
          annotation_out_ch.gff.join( annotation_out_ch.genome ),
          custom_db.mix( GET_NCBI_PROTEIN.out.proteins )
        )
        ch_custom_databases_annotations = CUSTOM_DATABASE.out.gff.groupTuple()
      }

      /*
          Eighth step -- Merge all annotations
      */
      MERGE_ANNOTATIONS( 
        annotation_out_ch.gff
          .join(kofamscan_output_ch,             remainder: true)
          .join(vfdb_output_ch,                  remainder: true)
          .join(victors_output_ch,               remainder: true)
          .join(amrfinder_output_ch,             remainder: true)
          .join(resfinder_gff_ch,                remainder: true)
          .join(rgi_output_ch,                   remainder: true)
          .join(iceberg_output_blastp_ch,        remainder: true)
          .join(phast_output_ch,                 remainder: true)
          .join(DIGIS.out.gff,                   remainder: true)
          .join(ch_custom_databases_annotations, remainder: true) 
      )

      /*
          Final step -- Create genome browser and reports' files
      */
      // Plot genomic islands
      DRAW_GIS( 
        MERGE_ANNOTATIONS.out.gff.join( ISLANDPATH.out.results ) 
      )

      // Convert GFF file to GBK file
      GFF2GBK( 
        MERGE_ANNOTATIONS.out.gff.join( annotation_out_ch.genome ) 
      )

      // Convert GFF file to sqldb
      CREATE_SQL(
        MERGE_ANNOTATIONS.out.gff
          .join( annotation_out_ch.genes     )
          .join( annotation_out_ch.proteins  )
          .join( annotation_out_ch.genome    )
          .join( DIGIS.out.gff_and_sequences )
      )

      JBROWSE(
        MERGE_ANNOTATIONS.out.gff
          .join( annotation_out_ch.genome                        )
          .join( annotation_out_ch.gff                           )
          .join( BARRNAP.out.gff                                 )
          .join( COMPUTE_GC.out.results                          )
          .join( resfinder_gff_ch,               remainder: true )
          .join( phigaro_output_bed_ch,          remainder: true )
          .join( ISLANDPATH.out.results,         remainder: true )
          .join( CALL_METHYLATION.out.bedgraph,  remainder: true )
          .join( CALL_METHYLATION.out.chr_sizes, remainder: true )
          .join( phispy_output_ch,               remainder: true )
          .join( MERGE_ANNOTATIONS.out.digis_gff                 )
          .join( antismash_output_ch,            remainder: true )
          .join( MERGE_ANNOTATIONS.out.customdb_gff.groupTuple(), remainder: true )
      )

      // Render reports
      if (params.custom_db || params.ncbi_proteins) {
        CUSTOM_DATABASE_REPORT( 
          CUSTOM_DATABASE.out.summary.join( MERGE_ANNOTATIONS.out.gff, remainder:true ) 
        )
      }
      REPORT(
        annotation_out_ch[8]
          .join( MERGE_ANNOTATIONS.out.gff                        )
          .join( BARRNAP.out.gff                                  )
          .join( MLST.out.results                                 )
          .join( kegg_decoder_svg_ch,             remainder: true )
          .join( REFSEQ_MASHER.out.results                        )
          .join( amrfinder_output_ch,             remainder: true )
          .join( rgi_output_ch,                   remainder: true )
          .join( rgi_output_parsed_ch,            remainder: true )
          .join( rgi_heatmap_ch,                  remainder: true )
          .join( argminer_output_ch,              remainder: true )
          .join( resfinder_output_tab_ch,         remainder: true )
          .join( resfinder_output_pointfinder_ch, remainder: true )
          .join( resfinder_phenotable_ch,         remainder: true )
          .join( vfdb_output_ch,                  remainder: true )
          .join( victors_output_ch,               remainder: true )
          .join( phigaro_output_tsv_ch,           remainder: true )
          .join( phispy_output_ch,                remainder: true )
          .join( iceberg_output_blastp_ch,        remainder: true )
          .join( iceberg_output_blastn_ch,        remainder: true )
          .join( plasmidfinder_output_ch,         remainder: true )
          .join( platon_output_ch,                remainder: true )
          .join( DRAW_GIS.out.example,            remainder: true )
          .join( phast_output_ch,                 remainder: true )
          .join( MERGE_ANNOTATIONS.out.digis_gff                  )
      )

}
