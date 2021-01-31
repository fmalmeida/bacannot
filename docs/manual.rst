.. _manual:

Manual
======

.. code-block:: bash

  # Get help in the command line
  nextflow run fmalmeida/bacannot --help

Parameters description
^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--use_tower``
     - N
     - False
     - Triggers the pipeline to be launched via nextflow tower

   * - ``--tower_token``
     - Y (if ``--use_tower``)
     - NA
     - Your nextflow tower token. Used to launch the pipeline in your nextflow tower account

   * - ``--outdir``
     - Y
     - output
     - Name of directory to store output values

   * - ``--threads``
     - N
     - 2
     - Number of threads to use

   * - ``--genome``
     - Y (if raw reads are not used)
     - NA
     - Genome(s) to be annotated in FASTA file. Mutually exclusively with the use of raw reads.

   * - ``--sreads_single``
     - N (Y if assembled genome is not used)
     - NA
     - Path to short unpaired reads.

   * - ``--sreads_paired``
     - N (Y if assembled genome is not used)
     - NA
     - Path to short paired reads

   * - ``--lreads``
     - N (Y if assembled genome is not used)
     - NA
     - Path to longreads (ONT or Pacbio)

   * - ``--lreads_type``
     - N (Y if longreads are used)
     - NA
     - Longreads is used? If so, from which tech it is? Options: [ 'nanopore' or 'pacbio' ]

   * - ``--bedtools_merge_distance``
     - N
     - 0
     - Minimum number of required overlapping bases to merge genes

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
     - Identity (%) threshold to be used when annotating AMR genes with ARGminer and AMRFinderPlus. CARD-RGI is not affected.

   * - ``--blast_resistance_mincov``
     - N
     - 90
     - Coverage (%) threshold to be used when annotating AMR genes with ARGminer and AMRFinderPlus. CARD-RGI is not affected.

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

   * - ``--resfinder_species``
     - N
     - NA
     - Activate the resfinder annotation process using the give species panel. Check them out in `their page <https://cge.cbs.dtu.dk/services/ResFinder/>`_.

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

   * - ``--nanopolish_fast5_dir``
     - N
     - NA
     - Path to directory containing fast5 files to be used to call methylation. If null, the analysis will be skipped

   * - ``--nanopolish_fastq_reads``
     - N
     - NA
     - Path to fastq reads (related to fast5 files) that will be used to call methylation. If null, the analysis will be skipped


All this parameters are configurable through a configuration file. We encourage users to use the configuration
file since it will keep your execution cleaner and more readable. See a :ref:`config` example.

Examples
""""""""

For a better understanding of the usage we provided a feel examples. See :ref:`examples`
