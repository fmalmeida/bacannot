.. _config:

Configuration File
""""""""""""""""""

To download a configuration file template users just need to run ``nextflow run fmalmeida/bacannot --get_config``

Using a config file your code is lot more clean and concise: ``nextflow run fmalmeida/bacannot -c [path-to-config]``

Default configuration:

.. code-block:: groovy

  /*
   * Configuration File to run fmalmeida/bacannot pipeline.
   */

  /*

                                                      Required Parameters.
                                              This parameters must always be set

  */
  params {
            /*
             * General parameters
             */

  // The input file formats are mutually exclusive. Users must choose between giving an
  // assembled genome or raw reads to the pipeline.
  // Input genome -- Always in FASTA format. Users can annotate more than one genome
  // at once by using glob patters, such as "*.fasta"
    genome = ''

  // Input raw reads -- Always in FASTQ format. When using raw reads, the fmalmeida/mpgap
  // is also required to be available. Understand that nextflow loads and uses the reads
  // in channels randomly so we can't guarantee that they will have the same order. Thus,
  // if using more than one NGS read type, you must give path to reads of one sample per
  // execution, otherwise the reads can end up mixed and confused. However, if using only
  // one NGS read type, only illumina paired ends for example, you can specify reads of
  // various samples in one run using glob patters such as "*{1,2}.fastq" since nextflow
  // will load all the pairs into the channel, all samples will be analysed, however, it
  // only works when using one NGS read type.
    sreads_single = '' // Path to unpaired illumina reads, if available for the sample
    sreads_paired = '' // Path to paired end illumina reads, if available for the sample
    lreads = '' // Path to longreads (ONT or Pacbio), if available for the sample
    lreads_type = '' // Longreads is used? If so, from which tech it is? Options: [ nanopore or pacbio ]

  // Main output folder name. More than one bacannot annotation can be redirected
  // to the same output parameter. It is good to keep related annotations together.
  // A subdirectory with the filename will be created inside this directory.
    outdir = 'output'

  // Number of threads to be used by each software
    threads = 2

  // Number of minimum overlapping base pairs required for merging
  // Negative values, such as -20, means the number of required overlapping bases for merging.
  // Positive values, such as 5, means the maximum distance accepted between features for merging.
  // By default (if Blank), this process is not executed. For execution the user needs to provide a value
    bedtools_merge_distance = ''

            /*
             * NF tower setup parameters
             */

    use_tower   = false
    tower_token = ''

            /*
             * Prokka optional parameters
             */

  // Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
    prokka_kingdom = ''

  // Translation table code. Must be set if the above is set.
  // Example: params.prokka_genetic.code = 11
    prokka_genetic.code = false

  // Use rnammer instead of Barrnap? False or True?
    prokka_use_rnammer = false

            /*
             * Resfinder optional parameter
             */

  // Species panel to be used when annotating with Resfinder. If blank,
  // it will not be executed. Must be identical (without the *) as written
  // in their webservice https://cge.cbs.dtu.dk/services/ResFinder/.
    resfinder_species = ''

            /*
             * Handling the execution of processes
             *
             * By default, all processes are executed. These
             * parameters tells wheter NOT to run a process.
             *
             * Which means: false will allow its execution
             * while true will create a barrier and skip a process.

  */
  // (NOT RUN?) Plasmids annotation (controls PlasmidFinder execution)
    not_run_plasmid_search = false

  // (NOT RUN?) General Virulence annotation (controls VFDB and Victors scan)
    not_run_virulence_search = false

  // (NOT RUN?) Resistance annotation (controls AMRfinder and RGI)
    not_run_resistance_search = false

  // (NOT RUN?) ICE annotation (controls ICEberg annotation)
    not_run_iceberg_search = false

  // (NOT RUN?) prophage annotation (controls PHAST and Phigaro)
    not_run_prophage_search = false

  // (NOT RUN?) KO (KEGG Orthology) annotation
    not_run_kofamscan = false

            /*
             * Annotation thresholds to be used when scanning specific databases and features
             * Select a combination of thresholds that is meaningful for your data. Some of
             * the databases are protein-only, others are nucleotide only. We cannnot control
             * that and the databases will be scanned either if blastp or blastn using these
             * thresholds described here.
             */

  // Identity threshold for plasmid annotation
    plasmids_minid = 90

  // Coverage threshold for plasmid annotation
    plasmids_mincov = 60

  // Virulence genes identity threshold
    blast_virulence_minid = 90

  // Virulence genes coverage threshold
    blast_virulence_mincov = 80

  // AMR genes identity threshold
    blast_resistance_minid= 90

  // AMR genes coverage threshold
    blast_resistance_mincov = 80

  // MGEs (ICEs and Phages) identity threshold
    blast_MGEs_minid = 65

  // MGEs (ICEs and Phages) coverage threshold
    blast_MGEs_mincov = 65

            /*
             * Configure optional Methylation annotation with nanopolish
             * If left blank, it will not be executed. When both parameters are set
             * it will automatically execute nanopolish to call methylation
             *
             * For using these parameters, the pipeline must be used with one sample at a time
             * since we can't guaratee the order the files are picked by nextflow.
             */

    nanopolish_fast5_dir = ''   // Path to directory containing FAST5 files
    nanopolish_fastq_reads = '' // Path to fastq files (file related to FAST5 files above)

  }

  /*
                                          Configuration of Nextflow Scopes
   */

  //Trace Report
  trace {
      enabled = false
      file = "${params.outdir}" + "/annotation_pipeline_trace.txt"
      fields = 'task_id,name,status,exit,realtime,cpus,%cpu,memory,%mem,rss'
  }

  //Timeline Report
  timeline {
      enabled = false
      file = "${params.outdir}" + "/annotation_pipeline_timeline.html"
  }

  //Complete Report
  report {
      enabled = false
      file = "${params.outdir}" + "/annotation_pipeline_nextflow_report.html"
  }

  /*
                        Setting NF tower configurations
  */
  if (params.use_tower) {
  env.TOWER_ACCESS_TOKEN = params.tower_token
  tower {
      accessToken = params.tower_token
      enabled = params.use_tower
  }
  }

  /*
                  Configuration of Docker images usage
                  DO NOT change any of those
  */

  // Docker permissions
  docker {
    enabled = true
    runOptions = '-u $(id -u):root'
  }

  // Queue limit
  executor.$local.queueSize = 1

  // specific images
  process {
      withLabel: 'main' {
          container = 'fmalmeida/bacannot:latest'
      }

      withLabel: 'renv' {
          container = 'fmalmeida/bacannot:renv'
      }

      withLabel: 'jbrowse' {
          container = 'fmalmeida/bacannot:jbrowse'
      }

      withLabel: 'kofam' {
          container = 'fmalmeida/bacannot:kofamscan'
      }

      withLabel: 'assembly' {
          container = 'fmalmeida/mpgap'
      }
  }
