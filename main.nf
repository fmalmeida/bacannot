#!/usr/bin/env nextflow
nextflow.preview.dsl=2

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

                                  General Parameters - Mandatory

    --outdir <string>                              Output directory name

    --threads <int>                                Number of threads to use

    --genome <string>                              User has as input only one genome. Set path to the
                                                   genome in FASTA file.

    --genome_fofn <string>                         CSV file (two values) containing the genome FASTA file
                                                   and its respective output prefix. One genome per line.
                                                   FAST path in first field, prefix value in the second.

    --bedtools_merge_distance                      Minimum number of overlapping bases for gene merge
                                                   using bedtools merge. Negative values, such as -20, means
                                                   the number of required overlapping bases for merging.
                                                   Positive values, such as 5, means the maximum distance
                                                   accepted for merging. By default, this process is not executed.
                                                   For execution the user needs to provide a value


                                  Prokka complementary parameters

    --prokka_kingdom <string>                      Prokka annotation mode. Possibilities (default 'Bacteria'):
                                                   Archaea|Bacteria|Mitochondria|Viruses

    --prokka_genetic_code <int>                    Genetic Translation code. Must be set if kingdom is not
                                                   default (in blank).

    --prokka_use_rnammer                           Tells prokka wheter to use rnammer instead of barrnap.


                                  Diamond (blastx) search parameters

    --diamond_virulence_minid                      Min. identity % for virulence annotation

    --diamond_virulence_mincov                     Min. gene coverage for virulence annotation

    --diamond_MGEs_minid                           Min. identity % for ICEs and prophage annotation

    --diamond_MGEs_mincov                          Min. query coverage for ICEs and prophage annotation


                                  Configure Optional processes

    --not_run_virulence_search                     Tells wheter you do not want to execute virulence annotation

    --not_run_resistance_search                    Tells wheter you do not want to execute resistance annotation

    --not_run_iceberg_search                       Tells wheter you do not want to execute ICE annotation

    --not_run_prophage_search                      Tells wheter you do not want to execute prophage annotation

    --not_run_kofamscan                            Tells wheter you do not want to execute KO annotation with kofamscan


                                  Configure Optional Pangenome analysis with Roary

    --roary_reference_genomes <string>             Used to set path to reference genomes to be used in the pangenome
                                                   analysis with Roary. Whenever set, the pipeline will automatically
                                                   execute Roary pangenome analysis. Example: "path/reference/*.fasta"
                                                   They must be all in one directory and they must no be links. They
                                                   must be the hard file.


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
nextflow run fmalmeida/bacannot --threads 3 --outDir kp25X --genome Kp31_BC08.contigs.fasta --bedtools_merge_distance -20 --prokka_center UNB --diamond_virulence_minid 90 --diamond_virulence_mincov 80 --diamond_MGEs_minid 70 --diamond_MGEs_mincov 60 --diamond_minimum_alignment_length 200 --virulence_search --vfdb_search --victors_search --resistance_search --ice_search --prophage_search --execute_kofamscan --nanopolish_fast5_dir fast5_pass --nanopolish_fastq_reads Kp31_BC08.fastq


""".stripIndent()
}

/*
 * Check if user needs help
 */

params.help = false
 // Show help emssage
 if (params.help){
   helpMessage()
   //file('work').deleteDir()
   //file('.nextflow').deleteDir()
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
params.genome = ''
params.genome_fofn = ''
// Prokka parameters
params.prokka_kingdom = ''
params.prokka_genetic_code = false
params.prokka_use_rnammer = false
// Blast parameters
params.diamond_virulence_minid = 90
params.diamond_virulence_mincov = 90
params.diamond_MGEs_minid = 65
params.diamond_MGEs_mincov = 65
params.diamond_minimum_alignment_length = 200
params.not_run_virulence_search = false
params.not_run_resistance_search = false
params.not_run_iceberg_search = false
params.not_run_prophage_search = false
params.not_run_kofamscan = false
params.roary_reference_genomes = false

/*
 * Define log message
 */

log.info "=============================================================="
log.info " Docker-based, fmalmeida/bacannot, Genome Annotation Pipeline "
log.info "=============================================================="
def summary = [:]
summary['Output dir']   = "${params.outdir}"
summary['Threads'] = params.threads
if (params.not_run_virulence_search == false) {
summary['Blast % ID - Virulence Genes'] = params.diamond_virulence_minid
summary['Blast query coverage - Virulence Genes'] = params.diamond_virulence_mincov
}
if (params.not_run_iceberg_search == false | params.not_run_prophage_search == false) {
summary['Blast % ID - ICEs or Phages'] = params.diamond_MGEs_minid
summary['Blast query coverage - ICEs or Phages'] = params.diamond_MGEs_mincov
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
  threads: params.threads, genome: params.genome, genome_fofn: params.genome_fofn)

// Virulence annotation with VFDB
include { vfdb } from './modules/virulence_scan_vfdb.nf' params(outdir: params.outdir,
  threads: params.threads, diamond_virulence_minid: params.diamond_virulence_minid,
  diamond_virulence_mincov: params.diamond_virulence_mincov)

// Virulence annotation with Victors
include { victors } from './modules/virulence_scan_victors.nf' params(outdir: params.outdir,
  threads: params.threads, diamond_virulence_minid: params.diamond_virulence_minid,
  diamond_virulence_mincov: params.diamond_virulence_mincov)

// Prophage annotation with PHAST
include { phast } from './modules/prophage_scan_phast.nf' params(outdir: params.outdir,
  threads: params.threads, diamond_MGEs_minid: params.diamond_MGEs_minid,
  diamond_MGEs_mincov: params.diamond_MGEs_mincov)

// Prophage annotation with PHIGARO
include { phigaro } from './modules/prophage_scan_phigaro.nf' params(outdir: params.outdir,
  threads: params.threads)

// ICE annotation with ICEberg db
include { iceberg } from './modules/ices_scan_iceberg.nf' params(outdir: params.outdir,
  threads: params.threads, diamond_MGEs_minid: params.diamond_MGEs_minid,
  diamond_MGEs_mincov: params.diamond_MGEs_mincov)

// Prophage annotation with PHIGARO
include { find_GIs } from './modules/IslandPath_DIMOB.nf' params(outdir: params.outdir)

/*
 * Define custom workflows
 */

// Only for a single genome
workflow bacannot_single_genome_nf {
  take:
    input_genome
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
      (params.not_run_kofamscan == false) ? kofamscan(prokka.out[4]) : ''
      (params.not_run_kofamscan == false) ? kegg_decoder(kofamscan.out[1]) : ''

      // Virulence search
      if (params.not_run_virulence_search == false) {
        // VFDB
        vfdb(prokka.out[5], prokka.out[3])
        // Victors db
        victors(prokka.out[4], prokka.out[3])
      }

      // Prophage search
      if (params.not_run_prophage_search == false) {
        // PHAST db
        phast(prokka.out[4], prokka.out[3])
        // Phigaro software
        phigaro(prokka.out[3])
      }

      // ICEs search
      if (params.not_run_iceberg_search == false) {
        // ICEberg db
        iceberg(prokka.out[5], prokka.out[4], prokka.out[3])
        // IslandPath software
        find_GIs(prokka.out[2])
      }
}

/*
 * Define main workflow
 */

workflow {

  // User has a single genome as input?
  if (params.genome) {
    bacannot_single_genome_nf(
      Channel.fromPath(params.genome)
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
