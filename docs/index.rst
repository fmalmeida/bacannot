.. Generic Archaeal and Bacterial Annotation (bacannot) pipeline documentation master file, created by
   sphinx-quickstart on Wed Nov 13 09:34:56 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Bacannot - A generic prokaryotic genome annotation pipeline
===========================================================

`bacannot <https://github.com/fmalmeida/bacannot>`_ is a pipeline developed with `Nextflow <https://www.nextflow.io/docs/latest/index.html>`_ and `Docker <https://www.docker.com/>`_. It was designed to provide an easy-to-use framework for performing a comprehensive annotation on prokaryotic genomes. Bacannot can annotate resistance genes, virulence factors, genomic islands, prophages, methylation and more. It wraps up the following tools and databases:

* [Prokka](https://github.com/tseemann/prokka)
* [barrnap](https://github.com/tseemann/barrnap)
* [mlst](https://github.com/tseemann/mlst)
* [Roary](https://github.com/sanger-pathogens/Roary)
* [Nanopolish](https://github.com/jts/nanopolish)
* [RGI](https://github.com/arpcard/rgi)
* [AMRFinderPlus](https://github.com/ncbi/amr/wiki)
* [diamond](https://github.com/bbuchfink/diamond)
* [VFDB](http://www.mgc.ac.cn/VFs/main.htm)
* [Victors](http://www.phidias.us/victors/)
* [ICEberg](http://db-mml.sjtu.edu.cn/ICEberg/)
* [PHAST](http://phast.wishartlab.com/)
* [Phigaro](https://github.com/bobeobibo/phigaro)
* [KofamScan](https://github.com/takaram/kofam_scan)
* [IslandPath-DIMOB](https://github.com/brinkmanlab/islandpath)
* [JBrowse](http://jbrowse.org/)

.. toctree::
   :hidden:

   installation
   manual

Support Contact
===============
Whenever a doubt arise feel free to contact me at <almeidafmarques@gmail.com>
