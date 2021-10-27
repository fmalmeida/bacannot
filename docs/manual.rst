.. _manual:

Manual
======

.. code-block:: bash

  # Get help in the command line
  nextflow run fmalmeida/bacannot --help

Parameters description
^^^^^^^^^^^^^^^^^^^^^^

Input files (single genome analysis)
""""""""""""""""""""""""""""""""""""

.. note::

  These parameters must only be used when annotating a single genome. If running the pipeline with more than 1 input
  genomes users must set them in the samplesheet YAML file as described in :ref:`samplesheet`.

.. note::
  
  Remember to always write input paths inside double quotes.

.. note::
  
  When using paired end reads it is required that input reads are set with the "{1,2}"" pattern. For example: "SRR6307304_{1,2}.fastq". This will properly load reads "SRR6307304_1.fastq" and "SRR6307304_2.fastq"

.. warning::
  
  When running hybrid assemblies or mixing short read types it is advised to **avoid not required REGEX** and write the full file path, using only the required REGEX for paired end reads when applicable. So that the pipeline does not load any different read that also matches the REGEX and avoid confusions with the inputs.

.. list-table::
   :widths: 20 10 20 30
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--genome``
     - Y (if raw reads are not used)
     - NA
     - Genome(s) to be annotated in FASTA file. Mutually exclusively with the use of raw reads.

   * - ``--prefix``
     - Y
     - out
     - This sets the prefix to be used when writing results

   * - ``--sreads_single``
     - N (Y if assembled genome is not used)
     - NA
     - Path to short unpaired reads. E.g. "SRR*.fastq.gz"

   * - ``--sreads_paired``
     - N (Y if assembled genome is not used)
     - NA
     - Path to short paired reads. E.g. "SRR6307304_{1,2}.fastq"

   * - ``--lreads``
     - N (Y if assembled genome is not used)
     - NA
     - Path to longreads (ONT or Pacbio)

   * - ``--lreads_type``
     - N (Y if longreads are used)
     - NA
     - Longreads are used? If so, from which technology it is? Options: [ 'nanopore' or 'pacbio' ]

Input files (multiple genome analysis)
""""""""""""""""""""""""""""""""""""""

.. list-table::
   :widths: 20 10 20 25
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--in_yaml``
     - Y
     - NA
     - Input samplesheet in YAML format. Used when analysis is to be performed with multiple genomes at once. It is incompatible with the parameters for single genome analysis.

Output directory
""""""""""""""""

.. list-table::
   :widths: 20 10 20 30
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--outdir``
     - Y
     - output
     - Name of directory to store output values. A sub-directory for each
       genome will be created inside this main directory.

Max job request
"""""""""""""""

.. list-table::
   :widths: 20 10 20 30
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--threads``
     - N
     - 2
     - Number of threads to use

   * - ``--parallel_jobs``
     - N
     - 1
     - Number of jobs to run in parallel. Each job can consume up to N threads (``--threads``)

Prokka annotation
"""""""""""""""""

.. list-table::
   :widths: 20 10 20 30
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

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
     - Tells Prokka whether to use rnammer instead of barrnap

Resfinder annotation
""""""""""""""""""""

.. note::

  This parameter must only be used when annotating a single genome. If running the pipeline in multi-sample mode,
  users must set it inside the samplesheet YAML file as described in :ref:`samplesheet`.

.. warning::

   Users must select one of the available Resfinder Species panels. They are listed at `their main page <https://cge.cbs.dtu.dk/services/ResFinder/>`_ and in `their repository page <https://bitbucket.org/genomicepidemiology/resfinder/src/master/#usage>`_. If your species is not available in Resfinder panels, you may use it with the "Other" panel (``--resfinder_species "Other"``).

.. list-table::
   :widths: 20 10 20 30
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--resfinder_species``
     - N
     - NA
     - Resfinder species panel. It activates the resfinder annotation process using the given species panel. Check the available species at `their page <https://cge.cbs.dtu.dk/services/ResFinder/>`_. If your species is not available in Resfinder panels, you may use it with the "Other" panel (``--resfinder_species "Other"``).

On/Off processes
""""""""""""""""

.. list-table::
   :widths: 20 10 20 30
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--skip_virulence_search``
     - N
     - False
     - Tells whether not to run virulence factors annotation. It skips both vfdb and victors annotation

   * - ``--skip_plasmid_search``
     - N
     - False
     - Tells whether not to run plasmid detection with Plasmidfinder

   * - ``--skip_resistance_search``
     - N
     - False
     - Tells whether not to run resistance genes annotation. It skips AMRFinderPlus and RGI annotation

   * - ``--skip_iceberg_search``
     - N
     - False
     - Tells whether not to run mobile genetic elements annotation with ICEberg

   * - ``--skip_prophage_search``
     - N
     - False
     - Tells whether not to run prophage annotation with PHAST and Phigaro

   * - ``--skip_kofamscan``
     - N
     - False
     - Tells whether not to run KEGG orthology (KO) annotation with KofamScan

   * - ``--skip_antismash``
     -  N
     - False
     - | Tells whether or not to run antiSMASH (secondary metabolite) annotation.
       | AntiSMASH is executed using only its core annotation modules in order to keep it fast

Custom nucl databases
"""""""""""""""""""""

.. list-table::
   :widths: 20 10 20 30
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--custom_db``
     - N
     - NA
     - Custom gene nucleotide databases to be used for additional annotations against the genome. See :ref:`custom-db`.

Annotation thresholds
"""""""""""""""""""""

.. list-table::
   :widths: 20 10 20 30
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--blast_virulence_minid``
     - N
     - 90
     - Identity (%) threshold to be used when annotating virulence factors from VFDB and Victors

   * - ``--blast_virulence_mincov``
     - N
     - 90
     - Coverage (%) threshold to be used when annotating virulence factors from VFDB and Victors

   * - ``--blast_resistance_minid``
     - N
     - 90
     - Identity (%) threshold to be used when annotating AMR genes with CARD-RGI, Resfinder, ARGminer and AMRFinderPlus.

   * - ``--blast_resistance_mincov``
     - N
     - 90
     - Coverage (%) threshold to be used when annotating AMR genes with Resfinder, ARGminer and AMRFinderPlus. CARD-RGI is not affected.

   * - ``--plasmids_minid``
     - N
     - 90
     - Identity (%) threshold to be used when detecting plasmids with Plasmidfinder

   * - ``--plasmids_mincov``
     - N
     - 60
     - Coverage (%) threshold to be used when detecting plasmids with Plasmidfinder

   * - ``--blast_MGEs_minid``
     - N
     - 85
     - Identity (%) threshold to be used when annotating prophages and mobile elements from PHAST and ICEberg databases

   * - ``--blast_MGEs_mincov``
     - N
     - 85
     - Coverage (%) threshold to be used when annotating prophages and mobile elements from PHAST and ICEberg databases

   * - ``--blast_custom_minid``
     - N
     - 0
     - Identity (%) threshold to be used when annotating with user's custom databases

   * - ``--blast_custom_mincov``
     - N
     - 0
     - Coverage (%) threshold to be used when annotating with user's custom databases

Methylation call
""""""""""""""""

.. note::

  This parameter must only be used when annotating a single genome. If running the pipeline with more than 1 input
  genomes users must set it in the samplesheet YAML file as described in :ref:`samplesheet`.

.. list-table::
   :widths: 20 10 20 30
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--nanopolish_fast5``
     - N
     - NA
     - Path to directory containing fast5 files to be used to call methylation. If null, the analysis will be skipped

   * - ``--nanopolish_fastq``
     - N
     - NA
     - Path to fastq reads (related to fast5 files) that will be used to call methylation. If null, the analysis will be skipped

Merge distance
""""""""""""""

.. list-table::
   :widths: 20 10 20 30
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--bedtools_merge_distance``
     - N
     - NA
     - Minimum number of required overlapping bases to merge genes. By default it is not executed.

All this parameters are configurable through a configuration file. We encourage users to use the configuration
file since it will keep your execution cleaner and more readable. See a :ref:`config` example.

Examples
^^^^^^^^

For a better understanding of the usage we provided a feel examples. See :ref:`examples`
