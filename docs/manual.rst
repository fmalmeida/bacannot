.. _manual:

Manual
======

.. code-block:: bash

  # Get help in the command line
  nextflow run fmalmeida/bacannot --help

Parameters description
^^^^^^^^^^^^^^^^^^^^^^

Input files
"""""""""""

Required
^^^^^^^^

To execute the annotation pipeline users **must** provide genomic data as either raw reads or assembled genomes as input. When raw reads are used, Unicycler and Flye assemblers are used to create, respectively, shortreads-only and hybrid assemblies, or longreads-only assemblies for the annotation process. Which means, the minimum required input files are:

* An assembled genome in FASTA format, **or**;
* Raw sequencing reads.

Optional
^^^^^^^^

The pipeline accepts as input two other input files types that are used to perform additional annotation processes, they are:

* path to a directory of FAST5

  * Then used together with nanopore reads it will call DNA methylation with Nanopolish.

* path to custom **nucleotide** databases as described in :ref:`custom-db`

  * These custom databases (``--custom_db``) will be used to perform additional annotation processes using BLASTn

.. note::

   Users must must carefully read the documentation in order to better understand the details of the pipeline workflow customization. Please read the :ref:`samplesheet manual page<samplesheet>` to better understand it.

Input samplesheet
"""""""""""""""""

.. list-table::
   :widths: 20 10 20 25
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--input``
     - Y
     - NA
     - Input samplesheet describing all the samples to be analysed.

Output directory
""""""""""""""""

.. list-table::
   :widths: 20 10 20 30
   :header-rows: 1

   * - Arguments
     - Required
     - Default value
     - Description

   * - ``--output``
     - Y
     - outdir
     - Name of directory to store output values. A sub-directory for each genome will be created inside this main directory.

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
     - NA
     - Number of jobs to run in parallel. Each job can consume up to N threads (``--threads``). If not given, let's nextflow automatically handle it.

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

  Sets a default value for input samples. If a sample has a different value given inside the samplesheet, the pipeline will use, for that sample, the value found inside the :ref:`samplesheet<samplesheet>`.

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
