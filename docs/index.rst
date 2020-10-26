.. Generic Archaeal and Bacterial Annotation (bacannot>`_ pipeline documentation master file, created by
   sphinx-quickstart on Wed Nov 13 09:34:56 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Bacannot - A generic genome annotation pipeline for prokaryotes
===============================================================

`bacannot <https://github.com/fmalmeida/bacannot>`_ is a pipeline developed with `Nextflow <https://www.nextflow.io/docs/latest/index.html>`_
and `Docker <https://www.docker.com/>`_. It was designed to provide an easy-to-use framework for performing a comprehensive annotation on
prokaryotic genomes. Bacannot can annotate resistance genes, virulence factors, genomic islands, prophages, methylation and more.
It wraps up the following tools and databases:

* `Prokka <https://github.com/tseemann/prokka>`_
* `barrnap <https://github.com/tseemann/barrnap>`_
* `mlst <https://github.com/tseemann/mlst>`_
* `Nanopolish <https://github.com/jts/nanopolish>`_
* `RGI <https://github.com/arpcard/rgi>`_
* `AMRFinderPlus <https://github.com/ncbi/amr/wiki>`_
* `diamond <https://github.com/bbuchfink/diamond>`_
* `VFDB <http://www.mgc.ac.cn/VFs/main.htm>`_
* `Victors <http://www.phidias.us/victors/>`_
* `ICEberg <http://db-mml.sjtu.edu.cn/ICEberg/>`_
* `PHAST <http://phast.wishartlab.com/>`_
* `Phigaro <https://github.com/bobeobibo/phigaro>`_
* `KofamScan <https://github.com/takaram/kofam_scan>`_
* `KEGGDecoder <https://github.com/bjtully/BioData/tree/master/KEGGDecoder>`_
* `IslandPath-DIMOB <https://github.com/brinkmanlab/islandpath>`_
* `Plasmidfinder <https://cge.cbs.dtu.dk/services/PlasmidFinder/>`_
* `JBrowse <http://jbrowse.org/>`_

.. toctree::
   :hidden:

   installation
   manual
   config
   examples

Support Contact
===============
Whenever a doubt arise feel free to contact me at almeidafmarques@gmail.com
