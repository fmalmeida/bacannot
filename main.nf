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
include { exampleMessage } from './nf_functions/examples.nf'
include { logMessage } from './nf_functions/log.nf'
include { write_csv } from './nf_functions/writeCSV.nf'
include { parse_csv } from './nf_functions/parseCSV.nf'
include { paramsCheck } from './nf_functions/paramsCheck.nf'
include { filter_ch } from './nf_functions/parseCSV.nf'


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

// Container manager
params.singularity = false
// General parameters
params.outdir = 'outdir'
params.threads = 2
params.bedtools_merge_distance = ''
// Input parameters
params.in_yaml = ''
params.prefix = ''
params.genome = ''
params.sreads_single = ''
params.sreads_paired = ''
params.lreads = ''
params.lreads_type = ''
// Prokka parameters
params.prokka_kingdom = ''
params.prokka_genetic_code = false
params.prokka_use_rnammer = false
// User custom db
params.custom_db = ''
params.blast_custom_minid = 0
params.blast_custom_mincov = 0
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
// Nanopolish
params.nanopolish_fast5 = ''
params.nanopolish_fastq = ''
// Workflow parameters
params.skip_plasmid_search = false
params.skip_virulence_search = false
params.skip_resistance_search = false
params.skip_iceberg_search = false
params.skip_prophage_search = false
params.skip_kofamscan = false

/*
 * Define log message
 */
logMessage()

/*
 * Include modules (Execution setup)
 */

// Unicycler assembly
include { unicycler_batch } from './modules/assembly/unicycler_batch.nf' params(outdir: params.outdir,
  threads: params.threads)

// Flye assembly
include { flye_batch } from './modules/assembly/flye_batch.nf' params(outdir: params.outdir,
  threads: params.threads)

// Unicycler assembly
include { unicycler } from './modules/assembly/unicycler.nf' params(outdir: params.outdir,
  threads: params.threads, prefix: params.prefix)

// Flye assembly
include { flye } from './modules/assembly/flye.nf' params(outdir: params.outdir,
  threads: params.threads, prefix: params.prefix, lreads_type: params.lreads_type)

/*
 * Define custom workflows
 */

// Parse samplesheet
include { parse_samplesheet } from './workflows/parse_samples.nf'

// Bacannot pipeline for multiple genomes
include { bacannot_nf } from './workflows/simple_workflow.nf'

// Bacannot pipeline for multiple genomes
include { bacannot_batch_nf } from './workflows/batch_workflow.nf'


/*
 * Define main workflow
 */

workflow {

  if (params.in_yaml) {

    parameter_yaml = new FileInputStream(new File(params.in_yaml))
    new Yaml().load(parameter_yaml).each { k, v -> params[k] = v }

    // Read YAML file
    parse_samplesheet(params.samplesheet)

    // Convert it to CSV for usability
    samples_ch = write_csv(parse_samplesheet.out)

    // Run unicycler when necessary
    unicycler_batch(filter_ch(samples_ch, "unicycler"))

    // Run flye when necessary
    flye_batch(filter_ch(samples_ch, "flye"))

    // Run annotation
    bacannot_batch_nf(filter_ch(samples_ch, "annotation").mix(flye_batch.out[1], unicycler_batch.out[1]),
                      (params.custom_db) ? Channel.fromPath( params.custom_db.split(',').collect{ it } ) : Channel.empty())

  } else {

    if (params.genome) {

      // User have an assembled genome
      bacannot_nf(Channel.fromPath(params.genome),
                 (params.nanopolish_fast5 && params.nanopolish_fastq) ? Channel.fromPath( params.nanopolish_fast5 )   : Channel.empty(),
                 (params.nanopolish_fast5 && params.nanopolish_fastq) ? Channel.fromPath( params.nanopolish_fastq )   : Channel.empty(),
                 (params.custom_db) ? Channel.fromPath( params.custom_db.split(',').collect{ it } ) : Channel.empty())

    } else if (params.sreads_single || params.sreads_paired) {

      // User have illumina reads (so it goes to unicycler)
      unicycler((params.sreads_paired) ? Channel.fromFilePairs( params.sreads_paired, flat: true, size: 2 ) : Channel.value(['', '', '']),
                (params.sreads_single) ? Channel.fromPath( params.sreads_single )                           : Channel.value(''),
                (params.lreads)        ? Channel.fromPath( params.lreads )                                  : Channel.value(''))
      bacannot_nf(unicycler.out[1],
                 (params.nanopolish_fast5 && params.nanopolish_fastq) ? Channel.fromPath( params.nanopolish_fast5 )   : Channel.empty(),
                 (params.nanopolish_fast5 && params.nanopolish_fastq) ? Channel.fromPath( params.nanopolish_fastq )   : Channel.empty(),
                 (params.custom_db) ? Channel.fromPath( params.custom_db.split(',').collect{ it } ) : Channel.empty())

    } else if ((params.lreads && params.lreads_type) && (!params.sreads_paired && !params.sreads_single)) {

      // User does not have illumina reads (so it goes to flye)
      flye(Channel.fromPath( params.lreads ))
      bacannot_nf(flye.out[1],
                 (params.nanopolish_fast5 && params.nanopolish_fastq) ? Channel.fromPath( params.nanopolish_fast5 )   : Channel.empty(),
                 (params.nanopolish_fast5 && params.nanopolish_fastq) ? Channel.fromPath( params.nanopolish_fastq )   : Channel.empty(),
                 (params.custom_db) ? Channel.fromPath( params.custom_db.split(',').collect{ it } ) : Channel.empty())
    }
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
