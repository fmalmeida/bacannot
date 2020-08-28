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

    --diamond_virulence_identity                   Min. identity % for virulence annotation

    --diamond_virulence_queryCoverage              Min. query coverage for virulence annotation

    --diamond_MGEs_identity                        Min. identity % for ICEs and prophage annotation

    --diamond_MGEs_queryCoverage                   Min. query coverage for ICEs and prophage annotation

    --diamond_minimum_alignment_length             Min. alignment length for diamond annotation


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
nextflow run fmalmeida/bacannot --threads 3 --outDir kp25X --genome Kp31_BC08.contigs.fasta --bedtools_merge_distance -20 --prokka_center UNB --diamond_virulence_identity 90 --diamond_virulence_queryCoverage 80 --diamond_MGEs_identity 70 --diamond_MGEs_queryCoverage 60 --diamond_minimum_alignment_length 200 --virulence_search --vfdb_search --victors_search --resistance_search --ice_search --prophage_search --execute_kofamscan --nanopolish_fast5_dir fast5_pass --nanopolish_fastq_reads Kp31_BC08.fastq


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
// Prokka parameters
params.prokka_kingdom = ''
params.prokka_genetic_code = false
params.prokka_use_rnammer = false
// Blast parameters
params.diamond_virulence_identity = 90
params.diamond_virulence_queryCoverage = 90
params.diamond_MGEs_identity = 65
params.diamond_MGEs_queryCoverage = 65
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
summary['Blast % ID - Virulence Genes'] = params.diamond_virulence_identity
summary['Blast query coverage - Virulence Genes'] = params.diamond_virulence_queryCoverage
}
if (params.not_run_iceberg_search == false | params.not_run_prophage_search == false) {
summary['Blast % ID - ICEs and Phages'] = params.diamond_MGEs_identity
summary['Blast query coverage - ICEs and Phages'] = params.diamond_MGEs_queryCoverage
}
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Configuration file'] = workflow.configFiles[0]
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 * Include modules (Execution setup)
 */

// Prokka annotation
include { prokka } from './modules/prokka.nf' params(outdir: params.outdir,
  prokka_kingdom: params.prokka_kingdom, prokka_genetic_code: params.prokka_genetic_code,
  prokka_use_rnammer: params.prokka_use_rnammer, threads: params.threads)

// MLST annotation
include { mlst } from './modules/mlst.nf' params(outdir: params.outdir)

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
      mlst(input_genome)
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
