.. _manual:

Manual
======

Overview
""""""""

.. image:: annotation_en.png

An overview of all annotation steps automatically taken by the pipeline.


Input
"""""

    * path to genome fasta file
    * path to referenge genomes fasta files
    * path to a directory of FAST5 files modified to contain basecall information
    * path to fastq reads

.. note::

   Users must **never** use hard or symbolic links. This will make nextflow fail.
   When setting the parameters, please **always** give full path to a hard file,
   not to a link. This will prevent file access fail.

Usage example
"""""""""""""

::

   nextflow run fmalmeida/bacannot [OPTIONS]

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--outDir``
     - Y
     - output
     - Name of directory to store output values

   * - ``--threads``
     - N
     - 2
     - Number of threads to use

   * - ``--genome``
     - Y
     - NA
     - Genome to be annotated in FASTA file

   * - ``--bedtools_merge_distance``
     - N
     - 0
     - Minimum number of required overlapping bases to merge genes

   * - ``--prokka_center``
     - Y
     - Centre
     - Your Institute acronym. It will be used by Prokka when renaming contigs.

   * - ``--prokka_kingdom``
     - N
     - Bacteria
     - Prokka annotation mode. Possibilities: Archaea|Bacteria

   * - ``--prokka_genetic_code``
     - N
     - 11
     - Genetic Translation code. Must be set if a different kingdom is customized.

   * - ``--prokka_use_rnammer``
     - N
     - False
     - Tells Prokka wheter to use rnammer instead of barrnap

   * - ``--prokka_genus``
     - N
     - NA
     - Set a specific prokka genus database to scan

   * - ``--diamond_virulence_identity``
     - N
     - 90
     - Identity (%) threshold to be used when annotating virulence factors from VFDB and Victors

   * - ``--diamond_virulence_queryCoverage``
     - N
     - 90
     - Coverage (%) threshold to be used when annotating virulence factors from VFDB and Victors

   * - ``--diamond_MGEs_identity``
     - N
     - 85
     - Identity (%) threshold to be used when annotating prophages and mobile elements from PHAST and ICEberg databases

   * - ``--diamond_MGEs_queryCoverage``
     - N
     - 85
     - Coverage (%) threshold to be used when annotating prophages and mobile elements from PHAST and ICEberg databases

   * - ``--diamond_minimum_alignment_length``
     - N
     - 200
     - Minimum alignment lenth to report a hit.

   * - ``--not_run_virulence_search``
     - N
     - False
     - Tells wheter not to run virulence factors annotation. It skips both vfdb and victors annotation

   * - ``--not_run_vfdb_search``
     - N
     - False
     - Tells wheter not to run virulence factors annotation with VFDB

   * - ``--not_run_victors_search``
     - N
     - False
     - Tells wheter not to run virulence factors annotation with Victors

   * - ``--not_run_resistance_search``
     - N
     - False
     - Tells wheter not to run resistance genes annotation. It skips AMRFinderPlus and RGI annotation

   * - ``--not_run_iceberg_search``
     - N
     - False
     - Tells wheter not to run mobile genetic elements annotation with ICEberg

   * - ``--not_run_prophage_search``
     - N
     - False
     - Tells wheter not to run prophage annotation with PHAST and Phigaro

   * - ``--not_run_kofamscan``
     - N
     - False
     - Tells wheter not to run KEGG orthology (KO) annotation with KofamScan

   * - ``--roary_reference_genomes``
     - N
     - NA
     - Path to reference genomes to be used in pangenome analysis. If null, the analysis will be skipped

   * - ``--nanopolish_fast5_dir``
     - N
     - NA
     - Path to directory containing fast5 files to be used to call methylation. If null, the analysis will be skipped

   * - ``--nanopolish_fastq_reads``
     - N
     - Path to fastq reads (related to fast5 files) that will be used to call methylation. If null, the analysis will be skipped





Configuration File
""""""""""""""""""

All this parameters are configurable through a configuration file. We encourage users to use the configuration file since it will keep your execution
cleaner and more readable.

To download a configuration file template users just need to run ``nextflow run fmalmeida/bacannot --get_config``

Using a config file your code is lot more clean and concise: ``nextflow run fmalmeida/bacannot -c [path-to-config]``

Exemplification of Config file:

.. code-block:: groovy

  /*

                                                      Required Parameters.
                                              This parameters must always be set

  */
  // Input file (Always in fasta)
  params.genome = ''

  // Name of output directory
  params.outDir = 'output'

  // Output files predix (it must never be written with white spaces)
  params.prefix = 'out'

  // Number of threads to be used
  params.threads = 2

  // Prokka will rename contig files into 'gnl|Centre|TAG_{1,2,3}'. Set your's institute acronym
  params.prokka_center = 'Centre'

  // Number of minimum overlapping base pairs required for merging
  // Negative values means required overlapping base pairs. Positive values means maximum distance accepted for merging.
  // Setting to false means using Bedtools default
  params.bedtools_merge_distance = false

  /*
   *
   *                                                 Prokka optional parameters
   *
   */
  // Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
  params.prokka_kingdom = ''

  // Translation table code. Must be set if the above is set.
  // Example: params.prokka_genetic.code = 11
  params.prokka_genetic.code = false

  // false or true to use rnammer instead of Barrnap
  params.prokka_use_rnammer = false

  // Set only if you want to search only a specific prokka genus database
  params.prokka_genus = ''

  /*
                                    DIAMOND parameters used to annotated the genome using the specific databases
                                    loaded in the docker image (VFDB, Victors, ICEberg and PHAST)
                                    Resistance and Virulence genes are searched with blastx while ICEs and Phages
                                    are searched with blastn.

  */
  // Virulence genes identity threshold
  params.diamond_virulence_identity = 90

  // Virulence genes coverage threshold
  params.diamond_virulence_queryCoverage = 90

  // MGEs (ICEs and Phages) identity threshold
  params.diamond_MGEs_identity = 85

  // MGEs (ICEs and Phages) coverage threshold
  params.diamond_MGEs_queryCoverage = 85

  // Minimum alignment length.
  params.diamond_minimum_alignment_length = 200

  /*

                                          Configure Optional Pangenome analysis with Roary
                                          Used to set path to reference genomes to be used in the pangenome
                                          analysis with Roary. Whenever set, the pipeline will automatically
                                          execute Roary pangenome analysis. Example: "path/reference/*.fasta"
                                          They must be all in one directory and they must no be links. They
                                          must be the hard file.

  */
  params.roary_reference_genomes = ''

  /*

                    Necessary files for calling methylation using nanopolish call-methylation algorithm.
                    This results will be readly plot in JBROWSE browser. Here we need Nanopore raw reads
                    and its fastq. This step is extremely time consuming. If you desire fast results it
                    is advised to skip this process and execute it later since it is not a difficult proccess.

                    To skip it one just need to left its variables blank.

  */
  params.fast5_dir = ''
  params.fastq_reads = ''

  /*

                                        Handling the execution of optional processes

                                        By default, all processes are executed. These
                                        parameters tells wheter NOT to run a process.

                                        Which means: false will allow its execution
                                        while true will create a barrier and skip a
                                        process.

  */
  // General Virulence annotation (this controls vfdb and victors together)
  params.not_run_virulence_search = false

  // Skip only VFDB annotation
  params.not_run_vfdb_search = false

  // Skip only Victors annotation
  params.not_run_victors_search = false

  // Skip Resistance annotation
  params.not_run_resistance_search = false

  // Skip ICE annotation
  params.not_run_iceberg_search = false

  // Skip prophage annotation
  params.not_run_prophage_search = false

  // Skip KO (KEGG Orthology) annotation
  params.not_run_kofamscan = false

Examples
""""""""

Check out the `examples <examples>`_ that we provide for a better understanding of the usage.
