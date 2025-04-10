/*
 * Include modules (Execution setup)
 */

include { UNICYCLER              } from '../modules/assembly/unicycler'
include { FLYE                   } from '../modules/assembly/flye'
include { REFSEQ_MASHER          } from '../modules/generic/mash'
include { SOURMASH_LCA           } from '../modules/generic/sourmash_lca'
include { SOURMASH_ALL           } from '../modules/generic/sourmash_all'
include { PROKKA                 } from '../modules/generic/prokka'
include { BAKTA                  } from '../modules/generic/bakta'
include { MLST                   } from '../modules/generic/mlst'
include { BARRNAP                } from '../modules/generic/barrnap'
include { COMPUTE_GC             } from '../modules/generic/compute_gc'
include { KOFAMSCAN              } from '../modules/KOs/kofamscan'
include { KEGG_DECODER           } from '../modules/KOs/kegg-decoder'
include { PLASMIDFINDER          } from '../modules/MGEs/plasmidfinder'
include { PLATON                 } from '../modules/MGEs/platon'
include { MOBSUITE               } from '../modules/MGEs/mob_suite'
include { VFDB                   } from '../modules/virulence/vfdb'
include { VICTORS                } from '../modules/virulence/victors'
include { PHAST                  } from '../modules/prophages/phast'
include { PHIGARO                } from '../modules/prophages/phigaro'
include { PHISPY                 } from '../modules/prophages/phispy'
include { ICEBERG                } from '../modules/MGEs/iceberg'
include { INTEGRON_FINDER        } from '../modules/MGEs/integron_finder'
include { INTEGRON_FINDER_2GFF   } from '../modules/MGEs/integron_finder_2gff'
include { ISLANDPATH             } from '../modules/MGEs/islandpath'
include { DRAW_GIS               } from '../modules/MGEs/draw_gis'
include { DIGIS                  } from '../modules/MGEs/digIS'
include { ARGMINER               } from '../modules/resistance/argminer'
include { RESFINDER              } from '../modules/resistance/resfinder'
include { AMRFINDER              } from '../modules/resistance/amrfinder'
include { CARD_RGI               } from '../modules/resistance/rgi_annotation'
include { CALL_METHYLATION       } from '../modules/generic/methylation'
include { CUSTOM_DATABASE        } from '../modules/generic/custom_database'
include { CUSTOM_DATABASE_REPORT } from '../modules/generic/custom_database_report'
include { GET_NCBI_PROTEIN       } from '../modules/generic/ncbi_protein'
include { GET_NCBI_GENOME        } from '../modules/generic/ncbi_genome'
include { MERGE_ANNOTATIONS      } from '../modules/generic/merge_annotations'
include { GFF2GBK                } from '../modules/generic/gff2gbk'
include { CREATE_SQL             } from '../modules/generic/gff2sql'
include { JBROWSE                } from '../modules/generic/jbrowse'
include { REPORT                 } from '../modules/generic/reports'
include { SEQUENCESERVER         } from '../modules/generic/sequenceserver'
include { ANTISMASH              } from '../modules/generic/antismash'
include { SUMMARY                } from '../modules/generic/summary'
include { MERGE_SUMMARIES        } from '../modules/generic/merge_summaries'
include { CIRCOS                 } from '../modules/generic/circos'

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
        kofamscan_all_ch    = KOFAMSCAN.out.all
        kofamscan_output_ch = KOFAMSCAN.out.results
        kegg_decoder_svg_ch = KEGG_DECODER.out.results
      } else {
        kofamscan_all_ch    = Channel.empty()
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
        plasmidfinder_all_ch    = PLASMIDFINDER.out.all
        plasmidfinder_output_ch = PLASMIDFINDER.out.results
        // platon
        PLATON( annotation_out_ch.genome, dbs_ch )
        platon_output_ch = PLATON.out.results
        platon_all_ch    = PLATON.out.all
        // mob suite
        MOBSUITE( annotation_out_ch.genome )
        mobsuite_output_ch = MOBSUITE.out.results
      } else {
        plasmidfinder_all_ch    = Channel.empty()
        plasmidfinder_output_ch = Channel.empty()
        platon_output_ch        = Channel.empty()
        platon_all_ch           = Channel.empty()
        mobsuite_output_ch      = Channel.empty()
      }

      // TODO: Maybe add in MGE optional?

      // IslandPath software
      ISLANDPATH( annotation_out_ch.gbk )

      // Integron_finder software
      if (!params.skip_integron_finder) {
        INTEGRON_FINDER( annotation_out_ch.genome )
        INTEGRON_FINDER_2GFF( INTEGRON_FINDER.out.gbk )
        ch_integron_finder_gff = INTEGRON_FINDER_2GFF.out.gff
      } else {
        ch_integron_finder_gff = Channel.empty()
      }

      // Virulence search
      if (params.skip_virulence_search == false) {     
        // VFDB
        VFDB( annotation_out_ch.genes, dbs_ch )
        vfdb_output_ch = VFDB.out.results
        vfdb_all_ch    = VFDB.out.all
        // Victors db
        VICTORS( annotation_out_ch.proteins, dbs_ch )
        victors_output_ch = VICTORS.out.results
        victors_all_ch    = VFDB.out.all
      } else {
        vfdb_all_ch       = Channel.empty()
        vfdb_output_ch    = Channel.empty()
        victors_output_ch = Channel.empty()
        victors_all_ch    = Channel.empty()
      }

      // Prophage search
      if (params.skip_prophage_search == false) {
        // PHAST db
        PHAST( annotation_out_ch.proteins, dbs_ch )
        phast_output_ch = PHAST.out.results
        phast_all_ch    = PHAST.out.all
        // Phigaro software
        PHIGARO( annotation_out_ch.genome, dbs_ch )
        phigaro_output_tsv_ch = PHIGARO.out.tsv
        phigaro_output_bed_ch = PHIGARO.out.bed
        phigaro_all_ch        = PHIGARO.out.all
        // PhiSpy
        PHISPY( annotation_out_ch.gbk )
        phispy_output_ch = PHISPY.out.results
        phispy_all_ch    = PHISPY.out.all
      } else {
        phast_all_ch          = Channel.empty()
        phast_output_ch       = Channel.empty()
        phigaro_all_ch        = Channel.empty()
        phigaro_output_tsv_ch = Channel.empty()
        phigaro_output_bed_ch = Channel.empty()
        phispy_all_ch         = Channel.empty()
        phispy_output_ch      = Channel.empty()
      }

      // ICEs search
      if (params.skip_iceberg_search == false) {
        // ICEberg db
        ICEBERG( annotation_out_ch.proteins, annotation_out_ch.genome, dbs_ch )
        iceberg_output_blastp_ch   = ICEBERG.out.results
        iceberg_output_blastn_ch   = ICEBERG.out.genome_summary
        iceberg_all_ch             = ICEBERG.out.all
      } else {
        iceberg_all_ch             = Channel.empty()
        iceberg_output_blastp_ch   = Channel.empty()
        iceberg_output_blastn_ch   = Channel.empty()
      }

      // AMR search
      if (params.skip_resistance_search == false) {
        // AMRFinderPlus
        AMRFINDER( annotation_out_ch.proteins, dbs_ch )
        amrfinder_output_ch = AMRFINDER.out.resistance_results
        amrfinder_all_ch    = AMRFINDER.out.all
        // CARD-RGI
        CARD_RGI( annotation_out_ch.proteins, dbs_ch )
        rgi_all_ch           = CARD_RGI.out.all
        rgi_output_ch        = CARD_RGI.out.raw_hits
        rgi_output_parsed_ch = CARD_RGI.out.parsed_hits
        rgi_heatmap_ch       = CARD_RGI.out.heatmap_png
        // ARGMiner
        ARGMINER( annotation_out_ch.proteins, dbs_ch )
        argminer_output_ch = ARGMINER.out.summary
        argminer_all_ch    = ARGMINER.out.all
        // Resfinder
        RESFINDER( annotation_out_ch.genome_with_species, dbs_ch )
        resfinder_all_ch                  = RESFINDER.out.all
        resfinder_output_tab_ch           = RESFINDER.out.results
        resfinder_output_pointfinder_ch   = RESFINDER.out.pointfinder_results
        resfinder_phenotable_ch           = RESFINDER.out.pheno_table
        resfinder_gff_ch                  = RESFINDER.out.gff
      } else {
        rgi_all_ch                        = Channel.empty()
        rgi_output_ch                     = Channel.empty()
        rgi_output_parsed_ch              = Channel.empty()
        rgi_heatmap_ch                    = Channel.empty()
        amrfinder_all_ch                  = Channel.empty()
        amrfinder_output_ch               = Channel.empty()
        argminer_all_ch                   = Channel.empty()
        argminer_output_ch                = Channel.empty()
        resfinder_all_ch                  = Channel.empty()
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

      //
      // BEGIN: sourmash-related modules
      //
      if (!params.skip_sourmash) {
        SOURMASH_LCA(
          dbs_ch,
          annotation_out_ch.genome,
          params.sourmash_scale,
          params.sourmash_kmer
        )

        // mashing against samples and close related genomes
        GET_NCBI_GENOME(
          REFSEQ_MASHER.out.results
          .map { it[1] }
          .splitCsv( sep: '\t', header: true )
          .map{ "${it.biosample}" }
          .unique()
        )

        SOURMASH_ALL(
          annotation_out_ch.genome
          .map{ it[1] }
          .mix( GET_NCBI_GENOME.out.genomes.collect() )
          .collect(),
          params.sourmash_scale,
          params.sourmash_kmer
        )

        ch_sourmash_plot = SOURMASH_ALL.out.plot.first()
      } else {
        ch_sourmash_plot = []
      }

      //
      // END: sourmash related modules
      //

      // IS identification
      DIGIS( annotation_out_ch.genome.join(annotation_out_ch.gbk) )

      // antiSMASH
      if (params.skip_antismash == false) {
        ANTISMASH( annotation_out_ch.gbk, dbs_ch )
        antismash_all_ch    = ANTISMASH.out.all
        antismash_output_ch = ANTISMASH.out.gff
      } else {
        antismash_all_ch    = Channel.empty()
        antismash_output_ch = Channel.empty()
      }

      // sequenceserver
      SEQUENCESERVER(
        annotation_out_ch.genome
          .join( annotation_out_ch.genes    )
          .join( annotation_out_ch.proteins )
      )

      // custom databases annotation
      if (params.custom_db || params.ncbi_proteins) {
        GET_NCBI_PROTEIN( ncbi_accs )
        CUSTOM_DATABASE(
          annotation_out_ch.gff.join( annotation_out_ch.genome ),
          custom_db.mix( GET_NCBI_PROTEIN.out.proteins )
        )
        ch_custom_annotations     = CUSTOM_DATABASE.out.gff.groupTuple()
        ch_custom_annotations_all = CUSTOM_DATABASE.out.all.groupTuple()
      } else {
        ch_custom_annotations     = Channel.empty()
        ch_custom_annotations_all = Channel.empty()
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
          .join(ch_custom_annotations,           remainder: true)
          .join(ch_integron_finder_gff,          remainder: true)
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
          .join( ch_integron_finder_gff,         remainder: true )
      )

      // Render reports
      if (params.custom_db || params.ncbi_proteins) {
        // parse GFFs
        custom_db_gffs_ch = 
        MERGE_ANNOTATIONS.out.customdb_gff
        .map{ id, file ->
            def db = file.baseName.replaceAll('^custom_database_', '')
            [ id, db, file ]
        }

        // report
        CUSTOM_DATABASE_REPORT( 
          CUSTOM_DATABASE.out.summary.join( custom_db_gffs_ch, by: [0, 1], remainder:true )
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
          .join( mobsuite_output_ch,              remainder: true )
          .join( DRAW_GIS.out.example,            remainder: true )
          .join( phast_output_ch,                 remainder: true )
          .join( MERGE_ANNOTATIONS.out.digis_gff                  )
          .join( ch_integron_finder_gff,          remainder: true ),
        ch_sourmash_plot
      )

      //
      // Generate annotation JSON summary
      //
      SUMMARY( 
        annotation_out_ch.all
        .join( MLST.out.all                , remainder: true )
        .join( BARRNAP.out.all             , remainder: true )
        .join( kofamscan_all_ch            , remainder: true )
        .join( plasmidfinder_all_ch        , remainder: true )
        .join( platon_all_ch               , remainder: true )
        .join( ISLANDPATH.out.results      , remainder: true )
        .join( vfdb_all_ch                 , remainder: true )
        .join( victors_all_ch              , remainder: true )
        .join( phast_all_ch                , remainder: true )
        .join( phigaro_all_ch              , remainder: true )
        .join( phispy_all_ch               , remainder: true )
        .join( iceberg_all_ch              , remainder: true )
        .join( amrfinder_all_ch            , remainder: true )
        .join( rgi_all_ch                  , remainder: true )
        .join( argminer_all_ch             , remainder: true )
        .join( resfinder_all_ch            , remainder: true )
        .join( CALL_METHYLATION.out.all    , remainder: true )
        .join( REFSEQ_MASHER.out.results   , remainder: true )
        .join( DIGIS.out.all               , remainder: true )
        .join( antismash_all_ch            , remainder: true )
        .join( MERGE_ANNOTATIONS.out.all   , remainder: true )
        .join( ch_integron_finder_gff      , remainder: true )
        .join( mobsuite_output_ch          , remainder: true )
      )
      MERGE_SUMMARIES(
        SUMMARY.out.summaries.map{ it[1] }.collect()
      )

      // Render circos plots
      if (!params.skip_circos) {
        circos_input_ch =
          annotation_out_ch.genome
          .join( annotation_out_ch.gff    , remainder: true )
          .join( MERGE_ANNOTATIONS.out.gff, remainder: true )
          .join( PHISPY.out.gff           , remainder: true )
          .map{
            it ->
              sample = it[0]
              it.remove(0)
              [ sample, it ]
          }
        CIRCOS( circos_input_ch )
      }

}
