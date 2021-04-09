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

              SINGLE GENOME ANALYSIS

      */

  // Prefix for writing genome assembly and annotatin resulting files
  // Preferentially the sample name
    prefix = 'out'

  // The input file formats are mutually exclusive. Users must choose between giving an
  // assembled genome or raw reads to the pipeline.
  // Input genome -- Always in FASTA format.
    genome = ''

  // Input raw reads -- Always in FASTQ format.
  // When using raw reads, the fmalmeida/mpgap is also required to be available.
    sreads_single = ''  // Path to unpaired illumina reads, if available for the sample
    sreads_paired = ''  // Path to paired end illumina reads, if available for the sample
    lreads = ''         // Path to longreads (ONT or Pacbio), if available for the sample
    lreads_type = ''    // Longreads is used? If so, from which tech it is? Options: [ nanopore or pacbio ]

  // Species panel to be used when annotating with Resfinder. If blank,
  // it will not be executed. Must be identical (without the *) as written
  // in their webservice https://cge.cbs.dtu.dk/services/ResFinder/.
  // E.g. 'Escherichia coli'; 'Klebsiella' ...
    resfinder_species = ''

  // Configure optional Methylation annotation with nanopolish
  // If left blank, it will not be executed. When both parameters are set
  // it will automatically execute nanopolish to call methylation

    nanopolish_fastq = ''   // Path to directory containing FAST5 files
    nanopolish_fast5 = ''   // Path to fastq files (file related to FAST5 files above)

      /*

              MULTIPLE GENOME ANALYSIS

      */

  // When analysing multiple genomes at once, all the parameters described above, must be, whenever
  // necessary and applicable to your data, set inside a samplesheet file in YAML format. We provide
  // an well-formated example of this YAML file at: https://github.com/fmalmeida/bacannot/blob/master/example_samplesheet.yaml
  //
  // Please read the example YAML samplesheet so you can understand how to properly fill it.
  //
  // It is also documented in the main manual: https://bacannot.readthedocs.io/en/latest/samplesheet.html
    in_yaml = ''

      /*

              GENERAL PARAMETERS -- FOR BOTH SINGLE AND MULTIPLE GENOME WORKFLOWS

       */

  // Main output folder name. More than one bacannot annotation can be redirected
  // to the same output parameter. It is good to keep related annotations together.
  // A subdirectory with the filename will be created inside this directory.
    outdir = 'output'

  // Number of threads to be used by each software
    threads = 2

  // Number of jobs to run in parallel. Be aware that each job (in parallel) can consume
  // N threads (set above). Be sure to carefully check your resources before augmenting
  // this parameter. For example: parallel_jobs = 2 + threads = 5 can consume up until 10
  // threads at once.
    parallel_jobs = 1

  // Number of minimum overlapping base pairs required for merging
  // Negative values, such as -20, means the number of required overlapping bases for merging.
  // Positive values, such as 5, means the maximum distance accepted between features for merging.
  // By default (if Blank), this process is not executed. For execution the user needs to provide a value
    bedtools_merge_distance = ''

            /*
             * Prokka optional parameters
             */

  // Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
    prokka_kingdom = ''

  // Translation table code. Must be set if the above is set.
  // Example: params.prokka_genetic.code = 11
    prokka_genetic_code = false

  // Use rnammer instead of Barrnap? False or True?
    prokka_use_rnammer = false

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
    skip_plasmid_search = false

  // (NOT RUN?) General Virulence annotation (controls VFDB and Victors scan)
    skip_virulence_search = false

  // (NOT RUN?) Resistance annotation (controls AMRfinder and RGI)
    skip_resistance_search = false

  // (NOT RUN?) ICE annotation (controls ICEberg annotation)
    skip_iceberg_search = false

  // (NOT RUN?) prophage annotation (controls PHAST and Phigaro)
    skip_prophage_search = false

  // (NOT RUN?) KO (KEGG Orthology) annotation
    skip_kofamscan = false

            /*
             * Custom databases can be used to annotate additional genes in the genome.
             * It runs a BLASTn alignment against the genome, therefore, the custom database
             * MUST be a nucleotide fasta of genes. More than one custom database can be given
             * separated by commas. Gene headers must be properly formated as described in the
             * documentation: https://bacannot.readthedocs.io/en/latest/custom-db.html
             */
  // Custom nucleotide fastas
    custom_db = ''

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

  // User's custom database identity threashold
    blast_custom_minid = 0

  // User's custom database coverage threashold
    blast_custom_mincov = 0

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
                  Setting up NF profiles
                  To use different profiles and executors
                  please read more at: https://www.nextflow.io/docs/latest/config.html#config-profiles
  */
  profiles {
    standard {
      // Executor
      process.executor = 'local'
      // QueueSize limit
      qs = (params.parallel_jobs) ? params.parallel_jobs : 1
      executor {
            name = 'local'
            queueSize = qs
      }
    }
  }
