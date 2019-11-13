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
    * path to a directory of FAST5 files modified to contain basecall information
    * path to fastq reads

.. note::

   Users must **never** use hard or symbolic links. This will make nextflow fail.
   When setting the parameters, please **always** give full path to a hard file,
   not to a link. This will prevent file access fail.

Usage example
"""""""""""""

::

   nextflow run fmalmeida/bacannot -c bacannot.config

.. list-table::
   :widths: 20 10 20 50
   :header-rows: 1

   * - Argument name(s)
     - Required
     - Default value
     - Description

   * -  <outDir>
     - Y
     - output
     - Name of directory to store output values

   * - ``--threads``
     - N
     - 2
     - Number of threads to use
