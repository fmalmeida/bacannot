#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import org.yaml.snakeyaml.Yaml

/*
 * Generic Pipeline for Prokariotic Genome Annotation
 */

/*
 * Include functions
 */
include { helpMessage } from './nf_functions/help.nf'
include { logMessage } from './nf_functions/log.nf'
include { write_csv } from './nf_functions/writeCSV.nf'
include { parse_csv } from './nf_functions/parseCSV.nf'
include { paramsCheck } from './nf_functions/paramsCheck.nf'


/*
 * Check parameters
 */
paramsCheck()
params.help = false
 // Show help emssage
 if (params.help){
   helpMessage()
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
params.output                  = 'outdir'
params.threads                 = 2
params.bedtools_merge_distance = ''
// Input parameters
params.input  = ''
// Prokka parameters
params.prokka_kingdom      = ''
params.prokka_genetic_code = false
params.prokka_use_rnammer  = false
// User custom db
params.custom_db           = ''
params.blast_custom_minid  = 0
params.blast_custom_mincov = 0
// Resfinder parameters
params.resfinder_species = ''
// Blast parameters
params.plasmids_minid          = 90
params.plasmids_mincov         = 60
params.blast_virulence_minid   = 90
params.blast_virulence_mincov  = 80
params.blast_resistance_minid  = 90
params.blast_resistance_mincov = 80
params.blast_MGEs_minid        = 65
params.blast_MGEs_mincov       = 65
// Workflow parameters
params.skip_plasmid_search    = false
params.skip_virulence_search  = false
params.skip_resistance_search = false
params.skip_iceberg_search    = false
params.skip_prophage_search   = false
params.skip_kofamscan         = false
params.skip_antismash         = false

/*
 * Define log message
 */
logMessage()

/*
 * Define custom workflows
 */

// Parse samplesheet
include { parse_samplesheet } from './workflows/parse_samples.nf'

// Bacannot pipeline for multiple genomes
include { BACANNOT } from './workflows/bacannot.nf'


/*
 * Define main workflow
 */

workflow {

  if (params.input) {

    // Load yaml
    parameter_yaml = file(params.input).readLines().join("\n")
    new Yaml().load(parameter_yaml).each { k, v -> params[k] = v }

    // Copy YAML samplesheet to output directory so user has a copy of it
    file(params.output).mkdir()
    file(params.input).copyTo(params.output + "/" + params.input)

    // Parse YAML file
    parse_samplesheet(params.samplesheet)

    // Convert it to CSV for usability
    samples_ch = write_csv(parse_samplesheet.out)

    // Run annotation
    BACANNOT(
      samples_ch,
      (params.custom_db) ? Channel.fromPath( params.custom_db.split(',').collect{ it } ) : Channel.empty()
    )

  } else {

    // Message to user
    println("""
    ERROR!
    A major error has occurred!
      ==> A samplesheet has not been provided. Please, provide a samplesheet to run the analysis. Online documentation is available at: https://bacannot.readthedocs.io/en/latest/
    Please, read the docs.
    Cheers.
    """)
  
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
