.. Generic Archaeal and Bacterial Annotation (bacannot)>`_

Bacannot - A generic genome annotation pipeline for prokaryotes
===============================================================

Designed to provide an easy-to-use framework for performing a comprehensive annotation on
prokaryotic genomes, `bacannot <https://github.com/fmalmeida/bacannot>`_ is developed with `Nextflow <https://www.nextflow.io/docs/latest/index.html>`_
and `Docker <https://www.docker.com/>`_. It can annotate resistance genes, virulence factors, genomic islands, prophages, methylation and more.

Its main steps are:

.. list-table::
   :widths: 60 40
   :header-rows: 1

   * - Analysis steps
     - Used software or databases

   * - Genome assembly (if raw reads are given)
     - `Flye <https://github.com/fenderglass/Flye>`_ and `Unicycler <https://github.com/rrwick/Unicycler>`_

   * - Generic annotation and gene prediction
     - `Prokka <https://github.com/tseemann/prokka>`_

   * - rRNA prediction
     - `barrnap <https://github.com/tseemann/barrnap>`_

   * - Classification within multi-locus sequence types (STs)
     - `mlst <https://github.com/tseemann/mlst>`_

   * - KEGG KO annotation and visualization
     - `KofamScan <https://github.com/takaram/kofam_scan>`_ and `KEGGDecoder <https://github.com/bjtully/BioData/tree/master/KEGGDecoder>`_

   * - Methylation annotation
     - `Nanopolish <https://github.com/jts/nanopolish>`_

   * - Annotation of antimicrobial (AMR) genes
     - `AMRFinderPlus <https://github.com/ncbi/amr/wiki>`_, `ARGminer <https://bench.cs.vt.edu/argminer>`_, `Resfinder <https://cge.cbs.dtu.dk/services/ResFinder/>`_ and `RGI <https://github.com/arpcard/rgi>`_

   * - Annotation of virulence genes
     - `Victors <http://www.phidias.us/victors/>`_ and `VFDB <http://www.mgc.ac.cn/VFs/main.htm>`_

   * - Prophage sequences and genes annotation
     - `PHASTER <http://phast.wishartlab.com/>`_, `Phigaro <https://github.com/bobeobibo/phigaro>`_ and `PhySpy <https://github.com/linsalrob/PhiSpy>`_

   * - Annotation of integrative and conjugative elements
     - `ICEberg <http://db-mml.sjtu.edu.cn/ICEberg/>`_

   * - *In silico* detection of plasmids
     - `Plasmidfinder <https://cge.cbs.dtu.dk/services/PlasmidFinder/>`_ and `Platon <https://github.com/oschwengers/platon>`_

   * - Prediction and visualization of genomic islands
     - `IslandPath-DIMOB <https://github.com/brinkmanlab/islandpath>`_ and `gff-toolbox <https://github.com/fmalmeida/gff-toolbox>`_

   * - Merge of annotation results
     - `bedtools <https://bedtools.readthedocs.io/en/latest/>`_

   * - Renderization of results in a Genome Browser
     - `JBrowse <http://jbrowse.org/>`_

   * - Renderization of automatic reports and shiny app for results interrogation
     - `R Markdown <https://rmarkdown.rstudio.com/>`_ and `Shiny <https://shiny.rstudio.com/>`_

.. toctree::
   :hidden:

   installation
   quickstart
   inputs
   outputs
   manual
   config
   custom-db
   samplesheet
   examples

Support Contact
===============
Whenever a doubt arise feel free to contact me at almeidafmarques@gmail.com
