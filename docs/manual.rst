.. _manual:

Manual
======

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

   * - ``--fast5_dir``
     - N
     - NA
     - Path to directory containing fast5 files to be used to call methylation. If null, the analysis will be skipped

   * - ``--fastq_reads``
     - N
     - NA
     - Path to fastq reads (related to fast5 files) that will be used to call methylation. If null, the analysis will be skipped


All this parameters are configurable through a configuration file. We encourage users to use the configuration
file since it will keep your execution cleaner and more readable. See a :ref:`config` example.

Examples
""""""""

For a better understanding of the usage we provided a feel examples. See :ref:`examples`
