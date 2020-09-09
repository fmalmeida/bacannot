.. _config:

Configuration File
""""""""""""""""""

To download a configuration file template users just need to run ``nextflow run fmalmeida/bacannot --get_config``

Using a config file your code is lot more clean and concise: ``nextflow run fmalmeida/bacannot -c [path-to-config]``

Default configuration:

.. code-block:: groovy

  /*

                                                      Required Parameters.
                                              This parameters must always be set

  */
  // Input genome -- Always in FASTA format
  params.genome = ''

  // Name of main output directory
  // A subdirectory with the filename will be created inside this directory
  params.outDir = 'output'

  // Number of threads to be used
  params.threads = 2

  // Number of minimum overlapping base pairs required for merging
  // Negative values, such as -20, means the number of required overlapping bases for merging.
  // Positive values, such as 5, means the maximum distance accepted between features for merging.
  // By default, this process is not executed. For execution the user needs to provide a value
  params.bedtools_merge_distance = ''

  /*

                                                    Prokka optional parameters

   */
  // Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
  params.prokka_kingdom = ''

  // Translation table code. Must be set if the above is set.
  // Example: params.prokka_genetic.code = 11
  params.prokka_genetic.code = false

  // false or true to use rnammer instead of Barrnap
  params.prokka_use_rnammer = false

  /*
                                    BLAST parameters used to annotated the genome using the specific databases
                                    loaded in the docker image (VFDB, Victors, ICEberg and PHAST) via blastx,
                                    blastp or blastn depending on the database characteristics.

  */
  // Virulence genes identity threshold
  params.blast_virulence_minid = 90

  // Virulence genes coverage threshold
  params.blast_virulence_mincov = 90

  // AMR genes identity threshold
  params.blast_resistance_minid= 90

  // AMR genes coverage threshold
  params.blast_resistance_mincov = 90

  // MGEs (ICEs and Phages) identity threshold
  params.blast_MGEs_minid = 85

  // MGEs (ICEs and Phages) coverage threshold
  params.blast_MGEs_mincov = 85


  /*
                                          Configure optional Methylation annotation with nanopolish
                                          If left blank, it will not be executed. When both parameters are set
                                          it will automatically execute nanopolish to call methylation
  */
  params.nanopolish_fast5_dir = ''   // Path to directory containing FAST5 files
  params.nanopolish_fastq_reads = '' // Path to fastq files (file related to FAST5 files above)

  /*

                                        Handling the execution of processes

                                        By default, all processes are executed. These
                                        parameters tells wheter NOT to run a process.

                                        Which means: false will allow its execution
                                        while true will create a barrier and skip a
                                        process.

  */
  // General Virulence annotation (controls VFDB scan)
  params.not_run_virulence_search = false

  // Skip Resistance annotation (controls AMRfinder and RGI)
  params.not_run_resistance_search = false

  // Skip ICE annotation (controls ICEberg annotation)
  params.not_run_iceberg_search = false

  // Skip prophage annotation (controls PHAST and Phigaro)
  params.not_run_prophage_search = false

  // Skip KO (KEGG Orthology) annotation
  params.not_run_kofamscan = false
