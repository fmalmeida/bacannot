#!/usr/bin/env nextflow
/*
========================================================================================
    fmalmeida/bacannot: A Generic Pipeline for Prokariotic Genome Annotation
========================================================================================
    Github : https://github.com/fmalmeida/bacannot
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2
import org.yaml.snakeyaml.Yaml

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOWS FOR PIPELINE
========================================================================================
*/

include { PARSE_SAMPLESHEET } from './workflows/parse_samples.nf'
include { BACANNOT          } from './workflows/bacannot.nf'
include { CREATE_DBS        } from './workflows/bacannot_dbs.nf'

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

workflow {

  if (params.get_dbs) {
    CREATE_DBS()
  } else {
    if (params.input) {

    // check if user gave path to bacannot databases
    if (!params.bacannot_db) {
      // Message to user
      exit("""
      ERROR!
      A major error has occurred!
        ==> User forgot to set path to databases with --bacannot_db. Online documentation is available at: https://bacannot.readthedocs.io/en/latest/
      Please, read the docs.
      Cheers.
      """)
    } else {
      bacannot_db = file(params.bacannot_db)
    }

    // Load yaml
    samplesheet_yaml = file(params.input)
    parameter_yaml = samplesheet_yaml.readLines().join("\n")
    new Yaml().load(parameter_yaml).each { k, v -> params[k] = v }

    // Copy YAML samplesheet to output directory so user has a copy of it
    file(params.output).mkdir()
    samplesheet_yaml.copyTo(params.output + "/" + "${samplesheet_yaml.getName()}")

    // Parse YAML file
    PARSE_SAMPLESHEET(params.samplesheet)

    // Run annotation
    BACANNOT(
      PARSE_SAMPLESHEET.out,
      bacannot_db,
      (params.custom_db) ? Channel.fromPath( params.custom_db.split(',').collect{ it } ) : Channel.empty(),
      (params.ncbi_proteins) ? Channel.fromPath( params.ncbi_proteins ) : Channel.empty()
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

}

/*
========================================================================================
    THE END
========================================================================================
*/
