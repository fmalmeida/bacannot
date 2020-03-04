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

   This pipeline contains a massive amount of configuration variables and its usage as CLI parameters would
   cause the command to be huge.
   Therefore, it is extremely recommended to use the nextflow.config configuration file in order to make
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
    --genome <string>                              Query Genome file
    --bedtools_merge_distance                      Minimum number of overlapping bases for gene merge
                                                   using bedtools merge (default: 0)

                            Prokka complementary parameters

    --prokka_center <string>                       Your institude acronym to be used by prokka when
                                                   renaming contigs.
    --prokka_kingdom <string>                      Prokka annotation mode. Possibilities (default 'Bacteria'):
                                                   Archaea|Bacteria|Mitochondria|Viruses
    --prokka_genetic_code <int>                    Genetic Translation code. Must be set if kingdom is not
                                                   default (in blank).
    --prokka_use_rnammer                           Tells prokka wheter to use rnammer instead of barrnap.
    --prokka_genus <string>                        Set only if you want to search only a specific genus database

                            Diamond (blastx) search parameters

    --diamond_virulence_identity                   Min. identity % for virulence annotation
    --diamond_virulence_queryCoverage              Min. query coverage for virulence annotation
    --diamond_MGEs_identity                        Min. identity % for ICEs and prophage annotation
    --diamond_MGEs_queryCoverage                   Min. query coverage for ICEs and prophage annotation
    --diamond_minimum_alignment_length             Min. alignment length for diamond annotation

                            Configure Optional processes

    --not_run_virulence_search                     Tells wheter you want or not to execute virulence annotation
    --not_run_vfdb_search                          Tells wheter you want or not to used VFDB database for virulence
                                                   annotation. It is useless if virulence_search is not true
    --not_run_victors_search                       Tells wheter you want or not to used victors database for virulence
                                                   annotation. It is useless if virulence_search is not true
    --not_run_resistance_search                    Tells wheter you want or not to execute resistance annotation
    --not_run_iceberg_search                       Tells wheter you want or not to execute ICE annotation
    --not_run_prophage_search                      Tells wheter you want or not to execute prophage annotation
    --not_run_kofamscan                            Tells wheter you want or not to execute KO annotation with kofamscan

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
  new File("bacannot.config") << new URL ("https://github.com/fmalmeida/bacannot/raw/master/nextflow.config").getText()
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

params.prefix = 'out'
params.outdir = 'outdir'
params.threads = 2
params.bedtools_merge_distance = 0
params.prokka_center = 'Centre'
params.prokka_kingdom = ''
params.prokka_genetic_code = false
params.prokka_use_rnammer = false
params.prokka_genus = ''
params.diamond_virulence_identity = 90
params.diamond_virulence_queryCoverage = 90
params.diamond_MGEs_identity = 65
params.diamond_MGEs_queryCoverage = 65
params.diamond_minimum_alignment_length = 200
params.not_run_virulence_search = false
params.not_run_vfdb_search = false
params.not_run_victors_search = false
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
summary['Input fasta']  = params.genome
summary['Output prefix']   = params.prefix
summary['Output dir']   = "${params.outdir}"
summary['Number of threads used'] = params.threads
summary['Blast % ID - Virulence Genes'] = params.diamond_virulence_identity
summary['Blast query coverage - Virulence Genes'] = params.diamond_virulence_queryCoverage
summary['Blast % ID - ICEs and Phages'] = params.diamond_MGEs_identity
summary['Blast query coverage - ICEs and Phages'] = params.diamond_MGEs_queryCoverage
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

// MLST analysis
include mlst from './modules/mlst.nf' params(outdir: params.outdir, prefix: params.prefix)

// Prokka annotation
include prokka from './modules/prokka.nf' params(outdir: params.outdir, prefix: params.prefix,
  prokka_kingdom: params.prokka_kingdom, prokka_genetic_code: params.prokka_genetic_code,
  prokka_use_rnammer: params.prokka_use_rnammer, prokka_genus: params.prokka_genus)

// Barrnap rRNA sequence prediction
include barrnap from './modules/barrnap.nf' params(outdir: params.outdir, prefix: params.prefix)

// Genome masking task
include masking_genome from './modules/genome_mask.nf' params(prefix: params.prefix)

// Compute GC content
include compute_gc from './modules/compute_gc.nf'

/*
 * Define custom workflows
 */

// Analysis for only one genome
workflow single_genome_nf {
  get:
    genome
    threads
  main:
    prokka(genome, threads)
    mlst(prokka.out[3])
    barrnap(prokka.out[3])
    masking_genome(prokka.out[3], prokka.out[1])
    compute_gc(prokka.out[3])
}

/*
 * Define main workflow
 */

workflow {
  // User gives only one genome
  if (params.genome) {
    single_genome_nf(Channel.fromPath(params.genome), params.threads)
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
