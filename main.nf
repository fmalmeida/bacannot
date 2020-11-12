#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Generic Pipeline for Prokariotic Genome Annotation
 */

/*
 * Define help message
 */

def helpMessage() {
   log.info """
   Usage:
   nextflow run fmalmeida/bacannot [--help] [ -c nextflow.config ] [OPTIONS] [-with-report] [-with-trace] [-with-timeline]
   Comments:

   This pipeline contains a massive amount of configuration variables and its usage as CLI parameters would cause the command
   to be huge. Therefore, it is extremely recommended to use the nextflow.config configuration file in order to make
   parameterization easier and more readable.

   Creating a configuration file:

   nextflow run fmalmeida/bacannot [--get_config]

   Show command line examples:

   nextflow run fmalmeida/bacannot --examples

   Execution Reports:

   nextflow run fmalmeida/bacannot [ -c nextflow.config ] -with-report
   nextflow run fmalmeida/bacannot [ -c nextflow.config ] -with-trace
   nextflow run fmalmeida/bacannot [ -c nextflow.config ] -with-timeline

   OBS: These reports can also be enabled through the configuration file.
   OPTIONS:

                                  General Parameters

    --outdir <string>                              Output directory name

    --threads <int>                                Number of threads to use

    --bedtools_merge_distance                      By default, this process is not executed. For execution
                                                   one needs to provide a value.Minimum number of overlapping
                                                   bases for gene merge using bedtools merge. Negative values,
                                                   such as -20, means the number of required overlapping bases
                                                   for merging. Positive values, such as 5, means the maximum
                                                   distance accepted between features for merging.

                                  Configuring Input options
                    (Users can give either a genome in FASTA file or raw reads)

    --genome <string>                              User has as input only one genome. Set path to the
                                                   genome in FASTA file.

    --sreads_paired                                Illumina paired end reads in fastq for assembly before annotation.

    --sreads_single                                Illumina unpaired reads in fastq for assembly before annotation.

    --lreads                                       Path to longreads in fastq assembly before annotation (ONT or Pacbio).

    --lreads_type                                  Tells the technology of the input longreads: [ nanopore or pacbio ].


                                  Prokka complementary parameters

    --prokka_kingdom <string>                      Prokka annotation mode. Possibilities (default 'Bacteria'):
                                                   Archaea|Bacteria|Mitochondria|Viruses

    --prokka_genetic_code <int>                    Genetic Translation code. Must be set if kingdom is not
                                                   default (in blank).

    --prokka_use_rnammer                           Tells prokka wheter to use rnammer instead of barrnap.


                                  Blast alignment parameters

    --blast_virulence_minid                        Min. identity % for virulence annotation. Default 90.

    --blast_virulence_mincov                       Min. gene coverage for virulence annotation. Default 80.

    --blast_resistance_minid                       Min. identity % for resistance annotation. Default 90.

    --blast_resistance_mincov                      Min. gene coverage for resistance annotation. Default 80.

    --blast_MGEs_minid                             Min. identity % for ICEs and prophage annotation. Default 65.

    --blast_MGEs_mincov                            Min. query coverage for ICEs and prophage annotation. Default 65.

    --plasmids_minid                               Min. identity % for plasmid detection. Default 90.

    --plasmids_mincov                              Min. query coverage for plasmid detection. Default 60.


                                  Configure resfinder optional parameter

    --resfinder_species                            It sets the species to be used for Resfinder annotation. If blank,
                                                   it will not be executed. Must be identical (without the *) as written
                                                   in their webservice https://cge.cbs.dtu.dk/services/ResFinder/.

                                  Configure Optional processes

    --not_run_virulence_search                     Tells wheter you do not want to execute virulence annotation

    --not_run_resistance_search                    Tells wheter you do not want to execute resistance annotation

    --not_run_iceberg_search                       Tells wheter you do not want to execute ICE annotation

    --not_run_prophage_search                      Tells wheter you do not want to execute prophage annotation

    --not_run_plasmid_search                       Tells wheter you do not want to execute plasmid detection

    --not_run_kofamscan                            Tells wheter you do not want to execute KO annotation with kofamscan


                            Configure optional Methylation annotation with nanopolish
                            If left blank, it will not be executed. And, with both parameters are set
                            it will automatically execute nanopolish to call methylation

    --nanopolish_fast5_dir <string>                Path to directory containing FAST5 files

    --nanopolish_fastq_reads <string>              Path to fastq files (file related to FAST5 files above)


""".stripIndent()
}

def exampleMessage() {
   log.info """
   Example Usages:
      Simple Klebsiella genome annotation using all pipeline's optional annotation processes
\$ nextflow run fmalmeida/bacannot --threads 3 --outdir kp25X --genome kp_ont.contigs.fasta --bedtools_merge_distance -20 --blast_virulence_minid 90 \
--blast_virulence_mincov 80 --blast_MGEs_minid 70 --blast_MGEs_mincov 60 --nanopolish_fast5_dir ./fast5_pass --nanopolish_fastq_reads ./kp_ont.fastq \
--resfinder_species "Klebsiella"


""".stripIndent()
}

/*
 * Check for errors
 */

// Input parameters

// Checking the use of lreads

// Prokka parameters
if (params.prokka_kingdom && !params.prokka_genetic_code) {
  log.info """
  ERROR!

  A minor error has occurred
    ==> User have set --prokka_kingdom but forget --prokka_genetic_code.

  These parameters must be used together. If you change prokka defaults kingdom parameter you must set the genetic code to be used for translation.

  If in doubt with these parameters let it blank, or get more information in Prokka's documentation.

  Cheers.
  """.stripIndent()

  exit 1
}

// Methylation parameters
if ((params.nanopolish_fast5_dir && !params.nanopolish_fastq_reads) || (!params.nanopolish_fast5_dir && params.nanopolish_fastq_reads)) {
  log.info """
  ERROR!

  A minor error has occurred
    ==> User have forget to set both --nanopolish_fast5_dir and --nanopolish_fastq_reads.

  These parameters must be used together. They are the necessary files to call methylations from ONT data with Nanopolish.

  Cheers.
  """.stripIndent()

  exit 1
}
/*
 * Check if user needs help
 */

params.help = false
 // Show help emssage
 if (params.help){
   helpMessage()
   exit 0
}

// CLI examples
params.examples = false
 // Show help emssage
 if (params.examples){
   exampleMessage()
   exit 0
}

/*
 * Does the user wants to download the configuration file?
 */

params.get_config = false
if (params.get_config) {
  new File("bacannot.config").write(new URL ("https://github.com/fmalmeida/bacannot/raw/master/nextflow.config").getText())
  println ""
  println "bacannot.config file saved in working directory"
  println "After configuration, run:"
  println "nextflow run fmalmeida/bacannot -c ./bacannot.config"
  println "Nice code!\n"
  exit 0
}

/*
 * Load general parameters and establish defaults
 */

// General parameters
params.outdir = 'outdir'
params.threads = 2
params.bedtools_merge_distance = ''
// Input parameters
params.genome = ''
params.sreads_single = ''
params.sreads_paired = ''
params.lreads = ''
params.lreads_type = ''
// Prokka parameters
params.prokka_kingdom = ''
params.prokka_genetic_code = false
params.prokka_use_rnammer = false
// Resfinder parameters
params.resfinder_species = ''
// Blast parameters
params.plasmids_minid = 90
params.plasmids_mincov = 60
params.blast_virulence_minid = 90
params.blast_virulence_mincov = 80
params.blast_resistance_minid = 90
params.blast_resistance_mincov = 80
params.blast_MGEs_minid = 65
params.blast_MGEs_mincov = 65
// Workflow parameters
params.not_run_plasmid_search = false
params.not_run_virulence_search = false
params.not_run_resistance_search = false
params.not_run_iceberg_search = false
params.not_run_prophage_search = false
params.not_run_kofamscan = false

/*
 * Define log message
 */

log.info "=============================================================="
log.info " Docker-based, fmalmeida/bacannot, Genome Annotation Pipeline "
log.info "=============================================================="
def summary = [:]
if (params.genome) { summary['Input genomes'] = params.genome }
summary['Output dir']   = "${params.outdir}"
summary['Threads'] = params.threads
if (params.not_run_virulence_search == false) {
summary['Blast % ID - Virulence Genes'] = params.blast_virulence_minid
summary['Blast query coverage - Virulence Genes'] = params.blast_virulence_mincov
}
if (params.not_run_resistance_search == false) {
summary['Blast % ID - AMR Genes'] = params.blast_resistance_minid
summary['Blast query coverage - AMR Genes'] = params.blast_resistance_mincov
}
if (params.not_run_iceberg_search == false | params.not_run_prophage_search == false) {
summary['Blast % ID - ICEs or Phages'] = params.blast_MGEs_minid
summary['Blast query coverage - ICEs or Phages'] = params.blast_MGEs_mincov
}
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Configuration file'] = workflow.configFiles[0]
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "=============================================================="

/*
 * Include modules (Execution setup)
 */

// Unicycler assembly
include { unicycler } from './modules/assembly/unicycler.nf' params(outdir: params.outdir,
  sreads_single: params.sreads_single, sreads_paired: params.sreads_paired,
  threads: params.threads, lreads: params.lreads)

// Flye assembly
include { flye } from './modules/assembly/flye.nf' params(outdir: params.outdir,
  threads: params.threads, lreads: params.lreads, lreads_type: params.lreads_type)

// Prokka annotation
include { prokka } from './modules/prokka.nf' params(outdir: params.outdir,
  prokka_kingdom: params.prokka_kingdom, prokka_genetic_code: params.prokka_genetic_code,
  prokka_use_rnammer: params.prokka_use_rnammer, threads: params.threads)

// MLST annotation
include { mlst } from './modules/mlst.nf' params(outdir: params.outdir)

// rRNA annotation
include { barrnap } from './modules/barrnap.nf' params(outdir: params.outdir)

// Calculate GC content
include { compute_gc } from './modules/compute_gc.nf'

// KOFAM annotation
include { kofamscan } from './modules/kofamscan.nf' params(outdir: params.outdir,
  threads: params.threads)

// KEGG decoder
include { kegg_decoder } from './modules/kegg-decoder.nf' params(outdir: params.outdir,
  threads: params.threads, genome: params.genome)

// Plasmid annotation with plasmidfinder
include { plasmidfinder } from './modules/plasmidfinder.nf' params(outdir: params.outdir,
  plasmids_minid: params.plasmids_minid, plasmids_mincov: params.plasmids_mincov)

// Virulence annotation with VFDB
include { vfdb } from './modules/virulence_scan_vfdb.nf' params(outdir: params.outdir,
  threads: params.threads, blast_virulence_minid: params.blast_virulence_minid,
  blast_virulence_mincov: params.blast_virulence_mincov)

// Virulence annotation with Victors
include { victors } from './modules/virulence_scan_victors.nf' params(outdir: params.outdir,
  threads: params.threads, blast_virulence_minid: params.blast_virulence_minid,
  blast_virulence_mincov: params.blast_virulence_mincov)

// Prophage annotation with PHAST
include { phast } from './modules/prophage_scan_phast.nf' params(outdir: params.outdir,
  threads: params.threads, blast_MGEs_minid: params.blast_MGEs_minid,
  blast_MGEs_mincov: params.blast_MGEs_mincov)

// Prophage annotation with PHIGARO
include { phigaro } from './modules/prophage_scan_phigaro.nf' params(outdir: params.outdir,
  threads: params.threads)

// ICE annotation with ICEberg db
include { iceberg } from './modules/ices_scan_iceberg.nf' params(outdir: params.outdir,
  threads: params.threads, blast_MGEs_minid: params.blast_MGEs_minid,
  blast_MGEs_mincov: params.blast_MGEs_mincov)

// Prophage annotation with PHIGARO
include { find_GIs } from './modules/IslandPath_DIMOB.nf' params(outdir: params.outdir)

// AMR annotation with ARGMiner
include { argminer } from './modules/resistance_scan_argminer.nf' params(outdir: params.outdir,
  threads: params.threads, blast_resistance_minid: params.blast_resistance_minid,
  blast_resistance_mincov: params.blast_resistance_mincov)

// AMR annotation with Resfinder
include { resfinder } from './modules/resistance_scan_resfinder.nf' params(outdir: params.outdir,
  threads: params.threads, blast_resistance_minid: params.blast_resistance_minid,
  blast_resistance_mincov: params.blast_resistance_mincov, resfinder_species: params.resfinder_species)

// AMR annotation with AMRFinderPlus
include { amrfinder } from './modules/amrfinder_scan.nf' params(outdir: params.outdir,
  threads: params.threads, blast_resistance_minid: params.blast_resistance_minid,
  blast_resistance_mincov: params.blast_resistance_mincov)

// AMR annotation with CARD-RGI
include { card_rgi } from './modules/rgi_annotation.nf' params(outdir: params.outdir,
  threads: params.threads, blast_resistance_minid: params.blast_resistance_minid)

// Methylation calling (Nanopolish)
include { call_methylation } from './modules/nanopolish_call_methylation.nf' params(outdir: params.outdir,
  threads: params.threads)

// Merging annotation in GFF
include { merge_annotations } from './modules/merge_annotations.nf' params(outdir: params.outdir)

// Convert GFF to GBK
include { gff2gbk } from './modules/gff2gbk.nf' params(outdir: params.outdir)

// Bedtools gff merge
include { gff_merge } from './modules/merge_gff.nf' params(outdir: params.outdir,
  bedtools_merge_distance: params.bedtools_merge_distance)

// JBrowse
include { jbrowse } from './modules/jbrowse.nf' params(outdir: params.outdir)

// MongoDB module
include { mongoDB } from './modules/create_mongoDB.nf' params(outdir: params.outdir)

// Output reports
include { report } from './modules/rmarkdown_reports.nf' params(outdir: params.outdir,
  blast_MGEs_mincov: params.blast_MGEs_mincov,
  blast_MGEs_minid: params.blast_MGEs_minid,
  blast_virulence_mincov: params.blast_virulence_mincov,
  blast_virulence_minid: params.blast_virulence_minid,
  blast_resistance_minid: params.blast_resistance_minid,
  blast_resistance_mincov: params.blast_resistance_mincov)

/*
 * Define custom workflows
 */

// Bacannot pipeline
workflow bacannot_nf {
  take:
    input_genome
    fast5_dir
    fast5_fastqs
  main:

      // First step -- Prokka annotation
      prokka(input_genome)

      // Second step -- MLST analysis
      mlst(prokka.out[3])

      // Third step -- rRNA annotation
      barrnap(prokka.out[3])

      // Fouth step -- calculate GC content for JBrowse
      compute_gc(prokka.out[3])

      // Fifth step -- run kofamscan
      if (params.not_run_kofamscan == false) {
        kofamscan(prokka.out[4])
        kegg_decoder(kofamscan.out[1])
        kofamscan_output = kofamscan.out[1]
      } else {
        kofamscan_output = Channel.empty()
      }

      // Sixth step -- MGE, Virulence and AMR annotations

      // Plasmid finder
      if (params.not_run_plasmid_search == false) {
        plasmidfinder(prokka.out[3])
        plasmidfinder_output = plasmidfinder.out[1]
      } else {
        plasmidfinder_output = Channel.empty()
      }

      // IslandPath software
      find_GIs(prokka.out[2])

      // Virulence search
      if (params.not_run_virulence_search == false) {
        // VFDB
        vfdb(prokka.out[5])
        vfdb_output = vfdb.out[1]
        // Victors db
        victors(prokka.out[4])
        victors_output = victors.out[1]
      } else {
        vfdb_output = Channel.empty()
        victors_output = Channel.empty()
      }

      // Prophage search
      if (params.not_run_prophage_search == false) {
        // PHAST db
        phast(prokka.out[4])
        phast_output = phast.out[1]
        // Phigaro software
        phigaro(prokka.out[3])
        phigaro_output = phigaro.out[1]
      } else {
        phast_output = Channel.empty()
        phigaro_output = Channel.empty()
      }

      // ICEs search
      if (params.not_run_iceberg_search == false) {
        // ICEberg db
        iceberg(prokka.out[4], prokka.out[3])
        iceberg_output = iceberg.out[1]
        iceberg_output_2 = iceberg.out[2]
      } else {
        iceberg_output = Channel.empty()
        iceberg_output_2 = Channel.empty()
      }

      // AMR search
      if (params.not_run_resistance_search == false) {
        // AMRFinderPlus
        amrfinder(prokka.out[4])
        amrfinder_output = amrfinder.out[0]
        // CARD-RGI
        card_rgi(prokka.out[4])
        rgi_output = card_rgi.out[3]
        rgi_output_1 = card_rgi.out[1]
        rgi_output_2 = card_rgi.out[2]
        // ARGMiner
        argminer(prokka.out[4])
        argminer_output = argminer.out[0]
        // Resfinder
        if (params.resfinder_species) {
          resfinder(prokka.out[3])
          resfinder_output_2 = resfinder.out[0]
          resfinder_output_1 = resfinder.out[1]
        } else {
          resfinder_output_1 = Channel.empty()
          resfinder_output_2 = Channel.empty()
        }
      } else {
        rgi_output = Channel.empty()
        rgi_output_1 = Channel.empty()
        rgi_output_2 = Channel.empty()
        amrfinder_output = Channel.empty()
        argminer_output = Channel.empty()
        resfinder_output_1 = Channel.empty()
        resfinder_output_2 = Channel.empty()
      }

      // Seventh step -- Methylation call
      if (params.nanopolish_fast5_dir && params.nanopolish_fastq_reads) {
        call_methylation(prokka.out[3], fast5_dir, fast5_fastqs)
        methylation_out_1 = call_methylation.out[2]
        methylation_out_2 = call_methylation.out[3]
      } else {
        methylation_out_1 = Channel.empty()
        methylation_out_2 = Channel.empty()
      }

      // Eighth step -- Merge all annotations with the same Prefix value in a single Channel
      annotations_files = prokka.out[3].join(prokka.out[1])
                                       .join(mlst.out[0])
                                       .join(barrnap.out[0])
                                       .join(compute_gc.out[0])
                                       .join(kofamscan_output, remainder: true)
                                       .join(vfdb_output, remainder: true)
                                       .join(victors_output, remainder: true)
                                       .join(amrfinder_output, remainder: true)
                                       .join(rgi_output, remainder: true)
                                       .join(iceberg_output, remainder: true)
                                       .join(phast_output, remainder: true)
                                       .join(phigaro_output, remainder: true)
                                       .join(find_GIs.out[0])

      // Contatenation of annotations in a single GFF file
      merge_annotations(annotations_files)

      // Convert GFF file to GBK file
      gff2gbk(merge_annotations.out[0].join(prokka.out[3]))

      // User wants to merge the final gff file?
      if (params.bedtools_merge_distance) {
        gff_merge(merge_annotations.out[0])
      }

      // Final step -- Create genome browser and reports

      // Grab inputs needed for JBrowse step
      jbrowse_input = merge_annotations.out[0].join(annotations_files, remainder: true)
                                              .join(methylation_out_1, remainder: true)
                                              .join(methylation_out_2, remainder: true)
      // Jbrowse Creation
      jbrowse(jbrowse_input)

      // Render reports
      report(jbrowse_input.join(rgi_output_1,         remainder: true)
                          .join(rgi_output_2,         remainder: true)
                          .join(argminer_output,      remainder: true)
                          .join(iceberg_output_2,     remainder: true)
                          .join(plasmidfinder_output, remainder: true)
                          .join(resfinder_output_1,   remainder: true)
                          .join(resfinder_output_2,   remainder: true))

}

/*
 * Define main workflow
 */

workflow {

  if (params.sreads_single || params.sreads_paired) {
    unicycler((params.sreads_paired) ? Channel.fromFilePairs( params.sreads_paired, flat: true, size: 2 ) : Channel.value(['', '', '']),
              (params.sreads_single) ? Channel.fromPath( params.sreads_single )                           : Channel.value(''),
              (params.lreads)        ? Channel.fromPath( params.lreads )                                  : Channel.value(''))
    bacannot_nf(unicycler.out[1],
               (params.nanopolish_fast5_dir && params.nanopolish_fastq_reads) ? Channel.fromPath( params.nanopolish_fast5_dir )   : Channel.empty(),
               (params.nanopolish_fast5_dir && params.nanopolish_fastq_reads) ? Channel.fromPath( params.nanopolish_fastq_reads ) : Channel.empty()
               )
  } else if ((params.lreads && params.lreads_type) && (!params.sreads_paired && !params.sreads_single)) {
    flye(Channel.fromPath( params.lreads ))
    bacannot_nf(flye.out[1],
               (params.nanopolish_fast5_dir && params.nanopolish_fastq_reads) ? Channel.fromPath( params.nanopolish_fast5_dir )   : Channel.empty(),
               (params.nanopolish_fast5_dir && params.nanopolish_fastq_reads) ? Channel.fromPath( params.nanopolish_fastq_reads ) : Channel.empty()
               )
  } else {
    bacannot_nf(Channel.fromPath(params.genome),
               (params.nanopolish_fast5_dir && params.nanopolish_fastq_reads) ? Channel.fromPath( params.nanopolish_fast5_dir )   : Channel.empty(),
               (params.nanopolish_fast5_dir && params.nanopolish_fastq_reads) ? Channel.fromPath( params.nanopolish_fastq_reads ) : Channel.empty()
               )
  }
}


// Completition message
workflow.onComplete {
    println ""
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Execution duration: $workflow.duration"
    println "Thank you for using fmalmeida/bacannot pipeline!"
}
