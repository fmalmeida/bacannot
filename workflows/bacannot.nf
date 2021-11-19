/*
 * Include modules (Execution setup)
 */

// Unicycler assembly
include { unicycler } from '../modules/assembly/unicycler.nf'

// Flye assembly
include { flye } from '../modules/assembly/flye.nf'

// Species identification
include { refseq_masher } from '../modules/generic/mash.nf'

// Prokka annotation
include { prokka } from '../modules/generic/prokka.nf'

// MLST annotation
include { mlst } from '../modules/generic/mlst.nf'

// rRNA annotation
include { barrnap } from '../modules/generic/barrnap.nf'

// Calculate GC content
include { compute_gc } from '../modules/generic/compute_gc.nf'

// KOFAM annotation
include { kofamscan } from '../modules/KOs/kofamscan.nf'

// KEGG decoder
include { kegg_decoder } from '../modules/KOs/kegg-decoder.nf'

// Plasmid annotation with plasmidfinder
include { plasmidfinder } from '../modules/MGEs/plasmidfinder.nf'

// Plasmid annotation with platon
include { platon } from '../modules/MGEs/platon.nf'

// Virulence annotation with VFDB
include { vfdb } from '../modules/virulence/vfdb.nf'

// Virulence annotation with Victors
include { victors } from '../modules/virulence/victors.nf'

// Prophage annotation with PHAST
include { phast } from '../modules/prophages/phast.nf'

// Prophage annotation with PHIGARO
include { phigaro } from '../modules/prophages/phigaro.nf'

// Prophage annotation with phispy
include { phispy } from '../modules/prophages/phispy.nf'

// ICE annotation with ICEberg db
include { iceberg } from '../modules/MGEs/iceberg.nf'

// Genomic Islands detection with Islandpath-DIMOB
include { find_GIs } from '../modules/MGEs/islandPath_DIMOB.nf'
include { draw_GIs } from '../modules/MGEs/draw_gis.nf'

// IS identification
include { digis } from '../modules/MGEs/digIS.nf'

// AMR annotation with ARGMiner
include { argminer } from '../modules/resistance/argminer.nf'

// AMR annotation with Resfinder
include { resfinder } from '../modules/resistance/resfinder.nf'

// AMR annotation with AMRFinderPlus
include { amrfinder } from '../modules/resistance/amrfinder.nf'

// AMR annotation with CARD-RGI
include { card_rgi } from '../modules/resistance/rgi_annotation.nf'

// Methylation calling (Nanopolish)
include { call_methylation } from '../modules/generic/methylation.nf'

// User's custom db annotation
include { custom_blast } from '../modules/generic/custom_blast.nf'
include { custom_blast_report } from '../modules/generic/custom_blast_report.nf'

// Merging annotation in GFF
include { merge_annotations } from '../modules/generic/merge_annotations.nf'

// Convert GFF to GBK
include { gff2gbk } from '../modules/generic/gff2gbk.nf'

// Convert GFF to SQL
include { create_sql } from '../modules/generic/gff2sql.nf'

// Bedtools gff merge
include { gff_merge } from '../modules/generic/merge_gff.nf'

// JBrowse
include { jbrowse } from '../modules/generic/jbrowse.nf'

// Output reports
include { report } from '../modules/generic/reports.nf'

// sequenceserver generation
include { sequenceserver } from '../modules/generic/sequenceserver.nf'

// antiSMASH
include { antismash } from '../modules/generic/antismash.nf'

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
        unicycler:  it[1] == 'unicycler'
        flye:       it[1] == 'flye'
        annotation: it[1] == "annotation"
      }.set { parsed_inputs }

      // Step 0 -- Run unicycler when necessary
      unicycler(parsed_inputs.unicycler)

      // Step 0 --  Run flye when necessary
      flye(parsed_inputs.flye)

      // First step -- Prokka or Bakta annotation
      if (params.use_bakta && params.bakta_db) {
        bakta(parsed_inputs.annotation.mix(flye.out[1], unicycler.out[1]))
        annotation_out_ch = bakta.out
      } else {
        prokka(parsed_inputs.annotation.mix(flye.out[1], unicycler.out[1]))
        annotation_out_ch = prokka.out
      }

      // Second step -- MLST analysis
      mlst(annotation_out_ch[3])

      // Third step -- rRNA annotation
      barrnap(annotation_out_ch[3])

      // Fouth step -- calculate GC content for JBrowse
      compute_gc(annotation_out_ch[3])

      // Fifth step -- run kofamscan
      if (params.skip_kofamscan == false) {
        kofamscan(annotation_out_ch[4])
        kegg_decoder(kofamscan.out[1])
        kofamscan_output = kofamscan.out[1]
        kegg_decoder_svg = kegg_decoder.out[1]
      } else {
        kofamscan_output = Channel.empty()
        kegg_decoder_svg = Channel.empty()
      }

      /*
          Sixth step -- MGE, Virulence and AMR annotations
      */

      // Plasmid finder
      if (params.skip_plasmid_search == false) {
        // plasmidfinder
        plasmidfinder(annotation_out_ch[3])
        plasmidfinder_output = plasmidfinder.out[1]
        // platon
        platon(annotation_out_ch[3])
        platon_output = platon.out[1]
      } else {
        plasmidfinder_output = Channel.empty()
        platon_output = Channel.empty()
      }

      // IslandPath software
      find_GIs(annotation_out_ch[2])

      // Virulence search
      if (params.skip_virulence_search == false) {
        // VFDB
        vfdb(annotation_out_ch[5])
        vfdb_output = vfdb.out[1]
        // Victors db
        victors(annotation_out_ch[4])
        victors_output = victors.out[1]
      } else {
        vfdb_output = Channel.empty()
        victors_output = Channel.empty()
      }

      // Prophage search
      if (params.skip_prophage_search == false) {
        // PHAST db
        phast(annotation_out_ch[4])
        phast_output = phast.out[1]
        // Phigaro software
        phigaro(annotation_out_ch[3])
        phigaro_output_1 = phigaro.out[0]
        phigaro_output_2 = phigaro.out[1]
        // PhiSpy
        phispy(annotation_out_ch[2])
        phispy_output = phispy.out[1]
      } else {
        phast_output = Channel.empty()
        phigaro_output_1 = Channel.empty()
        phigaro_output_2 = Channel.empty()
        phispy_output = Channel.empty()
      }

      // ICEs search
      if (params.skip_iceberg_search == false) {
        // ICEberg db
        iceberg(annotation_out_ch[4], annotation_out_ch[3])
        iceberg_output = iceberg.out[1]
        iceberg_output_2 = iceberg.out[2]
      } else {
        iceberg_output = Channel.empty()
        iceberg_output_2 = Channel.empty()
      }

      // AMR search
      if (params.skip_resistance_search == false) {
        // AMRFinderPlus
        amrfinder(annotation_out_ch[4])
        amrfinder_output = amrfinder.out[0]
        // CARD-RGI
        card_rgi(annotation_out_ch[4])
        rgi_output = card_rgi.out[2]
        rgi_output_parsed = card_rgi.out[1]
        rgi_heatmap = card_rgi.out[3]
        // ARGMiner
        argminer(annotation_out_ch[4])
        argminer_output = argminer.out[0]
        // Resfinder
        resfinder(annotation_out_ch[7])
        resfinder_output_1 = resfinder.out[0]
        resfinder_output_2 = resfinder.out[1]
        resfinder_phenotable = resfinder.out[2]
        resfinder_gff = resfinder.out[3]
      } else {
        rgi_output = Channel.empty()
        rgi_output_parsed = Channel.empty()
        rgi_heatmap = Channel.empty()
        amrfinder_output = Channel.empty()
        argminer_output = Channel.empty()
        resfinder_output_1 = Channel.empty()
        resfinder_output_2 = Channel.empty()
        resfinder_phenotable = Channel.empty()
        resfinder_gff = Channel.empty()
      }

      /*
          Seventh step -- Methylation call
      */
      call_methylation(annotation_out_ch[6])
      methylation_out_1 = call_methylation.out[2]
      methylation_out_2 = call_methylation.out[3]

      /*

          Additional steps created after main releases

       */

      // species identification
      refseq_masher(annotation_out_ch[3])

      // IS identification
      digis(annotation_out_ch[3].join(annotation_out_ch[2]))

      // antiSMASH
      if (params.skip_antismash == false) {
        antismash(annotation_out_ch[2])
        antismash_output = antismash.out[0]
      } else {
        antismash_output = Channel.empty()
      }

      // sequenceserver
      sequenceserver(annotation_out_ch[3].join(annotation_out_ch[5]).join(annotation_out_ch[4]))

      /*
          Eighth step -- Merge all annotations with the same Prefix value in a single Channel
      */
      annotations_files = annotation_out_ch[3].join(annotation_out_ch[1])
                                       .join(mlst.out[0])
                                       .join(barrnap.out[0])
                                       .join(compute_gc.out[0])
                                       .join(kofamscan_output, remainder: true)
                                       .join(vfdb_output,      remainder: true)
                                       .join(victors_output,   remainder: true)
                                       .join(amrfinder_output, remainder: true)
                                       .join(resfinder_gff,    remainder: true)
                                       .join(rgi_output,       remainder: true)
                                       .join(iceberg_output,   remainder: true)
                                       .join(phast_output,     remainder: true)
                                       .join(phigaro_output_2, remainder: true)
                                       .join(find_GIs.out[0],  remainder: true)

      // Contatenation of annotations in a single GFF file
      merge_annotations(annotations_files.join(digis.out[1],     remainder: true))

      // Plot genomic islands
      draw_GIs(merge_annotations.out[0].join(find_GIs.out[0]))

      // Convert GFF file to GBK file
      gff2gbk(merge_annotations.out[0].join(annotation_out_ch[3]))

      // Convert GFF file to sqldb
      create_sql(merge_annotations.out[0].join(annotation_out_ch[5])
                                         .join(annotation_out_ch[4])
                                         .join(annotation_out_ch[3])
                                         .join(digis.out[2]))

      // User wants to merge the final gff file?
      if (params.bedtools_merge_distance) {
        gff_merge(merge_annotations.out[0])
      }

      /*

          Nineth step -- Perform users custom annotation

      */
      if (params.custom_db) {
        custom_blast(merge_annotations.out[0].join(annotation_out_ch[3]), custom_db)
        custom_blast_report(custom_blast.out[0])
      }

      /*
          Final step -- Create genome browser and reports
      */

      // Grab inputs needed for JBrowse step
      jbrowse_input = merge_annotations.out[0].join(annotations_files, remainder: true)
                                              .join(methylation_out_1, remainder: true)
                                              .join(methylation_out_2, remainder: true)
                                              .join(phispy_output,     remainder: true)
                                              .join(merge_annotations.out[8], remainder: true) // parsed and changed digIS
                                              .join(antismash_output,  remainder: true)
      // Jbrowse Creation
      jbrowse(jbrowse_input)

      // Render reports
      report(jbrowse_input.join(rgi_output_parsed,    remainder: true)
                          .join(rgi_heatmap,          remainder: true)
                          .join(argminer_output,      remainder: true)
                          .join(iceberg_output_2,     remainder: true)
                          .join(plasmidfinder_output, remainder: true)
                          .join(resfinder_output_1,   remainder: true)
                          .join(resfinder_output_2,   remainder: true)
                          .join(resfinder_phenotable, remainder: true)
                          .join(draw_GIs.out[1],      remainder: true)
                          .join(phigaro_output_1,     remainder: true)
                          .join(platon_output,        remainder: true)
                          .join(annotation_out_ch[8],        remainder: true)
                          .join(kegg_decoder_svg,     remainder: true)
                          .join(refseq_masher.out[0], remainder: true))

}
