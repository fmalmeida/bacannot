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

params.outdir = 'outdir'
params.threads = 2
params.bedtools_merge_distance = ''
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
include mlst from './modules/mlst.nf' params(outdir: params.outdir)

// Prokka annotation
include prokka from './modules/prokka.nf' params(outdir: params.outdir,
  prokka_kingdom: params.prokka_kingdom, prokka_genetic_code: params.prokka_genetic_code,
  prokka_use_rnammer: params.prokka_use_rnammer, prokka_genus: params.prokka_genus,
  prokka_center: params.prokka_center, threads: params.threads)

// Barrnap rRNA sequence prediction
include barrnap from './modules/barrnap.nf' params(outdir: params.outdir)

// Genome masking task
// include masking_genome from './modules/genome_mask.nf'

// Compute GC content
include compute_gc from './modules/compute_gc.nf'

// Kofamscan analysis (Annotate KOs)
include kofamscan from './modules/kofamscan.nf' params(outdir: params.outdir, threads: params.threads)

// Virulence Annotation with VFDB
include vfdb from './modules/virulence_scan_vfdb.nf' params(outdir: params.outdir,
  diamond_virulence_queryCoverage: params.diamond_virulence_queryCoverage,
  diamond_virulence_identity: params.diamond_virulence_identity)

// Virulence Annotation with Victors
include victors from './modules/virulence_scan_victors.nf' params(outdir: params.outdir,
  diamond_virulence_queryCoverage: params.diamond_virulence_queryCoverage,
  diamond_virulence_identity: params.diamond_virulence_identity)

// Resistance annotation with AMRFINDERPLUS
include amrfinder from './modules/amrfinder_scan.nf' params(outdir: params.outdir)

// Resistance annotation with RGI
include rgi from './modules/rgi_annotation.nf' params(outdir: params.outdir, threads: params.threads)

// ICE annotation with ICEberg
include iceberg from './modules/ices_scan_iceberg.nf' params(outdir: params.outdir,
  diamond_MGEs_queryCoverage: params.diamond_MGEs_queryCoverage,
  diamond_MGEs_identity: params.diamond_MGEs_identity)

// PHAST annotation
include phast from './modules/prophage_scan_phast.nf' params(outdir: params.outdir,
  diamond_MGEs_queryCoverage: params.diamond_MGEs_queryCoverage,
  diamond_MGEs_identity: params.diamond_MGEs_identity)

// Prophage scan with phigaro
include phigaro from './modules/prophage_scan_phigaro.nf' params(outdir: params.outdir)

// GI prediction with IslandPath-DIMOB
include find_GIs from './modules/islandPath_DIMOB.nf' params(outdir: params.outdir,
  prokka_center: params.prokka_center)

// Merging annotation in GFF
include update_gff from './modules/updating_gff.nf' params(outdir: params.outdir)

// Convert GFF to GBK
include gff2gbk from './modules/gff2gbk.nf' params(outdir: params.outdir)

// Bedtools gff merge
include gff_merge from './modules/merge_gff.nf' params(outdir: params.outdir,
  bedtools_merge_distance: params.bedtools_merge_distance)

// Methylation calc (Nanopolish)
include call_methylation from './modules/nanopolish_call_methylation.nf' params(outdir: params.outdir)

/*
 * Define custom workflows
 */

// Analysis with methylation calling
workflow bacannot_nf {
  take:
    input_tuple
  main:
  // First we need to parse and separate the values in (main -> fasta and prefix) and (additional -> files for methylation annotation)
  input_tuple.multiMap { it ->
    base: [ it[0], it[1] ]
    add:  [ it[0], it[2], it[3] ]
  }.set { prokka_input }

  // Then we can start with the first process, Prokka annotation

  // Prokka Annotation
  prokka(prokka_input.base)

  // MLST analysis
  mlst(prokka.out[3])

  // rRNA prediction
  barrnap(prokka.out[3])

  // Genome masking and removal of nucleotide sequences from prokka gff
  // masking_genome(prokka.out[3], prokka.out[1], parsed_input.prefix)

  // Computing GC content
  compute_gc(prokka.out[3])

  // User wants kofamscan?
  if (params.not_run_kofamscan == false) {
    kofamscan(prokka.out[4])
    kofamscan_output = kofamscan.out[1]
  } else { kofamscan_output = Channel.empty().ifEmpty('EMPTY') }

  // User wants virulence search?
  if (params.not_run_virulence_search == false) {
    vfdb(prokka.out[5])
    vfdb_output = vfdb.out[0]
    victors(prokka.out[5])
    victors_output = victors.out[0]
  } else {
    vfdb_output    = Channel.empty().ifEmpty('EMPTY')
    victors_output = Channel.empty().ifEmpty('EMPTY')
  }

  // User wants resistance search?
  if (params.not_run_resistance_search == false) {
    amrfinder(prokka.out[4])
    amrfinder_output = amrfinder.out[0]
    rgi(prokka.out[4])
    rgi_output = rgi.out[2]
  } else {
    amrfinder_output = Channel.empty().ifEmpty('EMPTY')
    rgi_output       = Channel.empty().ifEmpty('EMPTY')
  }

  // User wants to use ICEberg db?
  if (params.not_run_iceberg_search == false) {
    iceberg(prokka.out[5])
    iceberg_output = iceberg.out[0]
  } else { iceberg_output = Channel.empty().ifEmpty('EMPTY') }

  // User wants prophage search?
  if (params.not_run_prophage_search == false) {
    phast(prokka.out[5])
    phast_output = phast.out[0]
    phigaro(prokka.out[3])
    phigaro_output = phigaro.out[1]
  } else {
    phast_output   = Channel.empty().ifEmpty('EMPTY')
    phigaro_output = Channel.empty().ifEmpty('EMPTY')
  }

  // Prediction of Genomic Islands
  find_GIs(prokka.out[2])

  // Merge all annotations with the same Prefix value in a single Channel
  annotations_files = prokka.out[3].join(prokka.out[1])
                                   .join(mlst.out[0])
                                   .join(barrnap.out[0])
                                   .join(compute_gc.out[0])
                                   .join(kofamscan_output)
                                   .join(vfdb_output)
                                   .join(victors_output)
                                   .join(amrfinder_output)
                                   .join(rgi_output)
                                   .join(iceberg_output)
                                   .join(phast_output)
                                   .join(phigaro_output)
                                   .join(find_GIs.out[0])

  // Contatenation of annotations in a single GFF file
  // update_gff(annotations_files)

  // Check this things
  annotations_files.view()

  // Convert GFF file to GBK file
  // gff2gbk(update_gff.out[0], prokka.out[3], parsed_input.prefix)

  // User wants to merge the final gff file?
  // if (params.bedtools_merge_distance != '') {
  //  gff_merge(update_gff.out[0], parsed_input.prefix)
  //}

  // Grab all input files with new genome/contig names for methylation call
  // methylation_input = prokka.out[3].join(prokka_input.add)

  // User wants to call methylation?
  // If not (based on input.fofn) the process will only create blank files
  // so the pipeline doesn't stop
  // call_methylation(methylation_input)

}

/*
 * Define main workflow
 */

workflow {
  // Channel for execution
  input = Channel.empty()

  // Read genomes, line by line
  fofn   = file(params.genome_fofn)
  lines  = fofn.readLines()
  for( line in lines ) {

    // User do not want methylation call

    if (line.split(',').size() == 2) {
      (prefix, fasta) = line.split(',');
      input << [prefix, file(fasta)]
    }

    // User wants methylation call

    else if (line.split(',').size() == 4) {
      // Execute the workflow for each value pair
      (prefix, fasta, fastq_reads, fast5_dir) = line.split(',');
      input << [prefix, file(fasta), file(fastq_reads), file(fast5_dir)]
    }
  }

  // Run
  bacannot_nf(input)
}


// Completition message
workflow.onComplete {
    println ""
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Execution duration: $workflow.duration"
    println "Thank you for using fmalmeida/bacannot pipeline!"
}
