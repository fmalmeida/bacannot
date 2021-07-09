.. _inputs:

Input files
===========

Required
^^^^^^^^

To execute the annotation pipeline users **must** provide genomic data as either raw reads or assembled genomes as input. When raw reads are used, Unicycler and Flye
assemblers are used to create, respectively, shortreads-only and hybrid assemblies, or longreads-only assemblies for the annotation process. Which means, the minimum
required input files are:

* An assembled genome in FASTA format, **or**;
* Raw sequencing reads.

.. note::

  Users can analyse more than one genome at once by properly configuring a samplesheet
  and using it with the ``--in_yaml`` parameter. See :ref:`samplesheet`.

Optional
^^^^^^^^

The pipeline accepts as input two other input files types that are used to perform additional annotation processes, they are:

* path to a directory of FAST5 and path to ONT fastq

  * This data will be used for the methylation calling process

* path to custom **nucleotide** databases as described in :ref:`custom-db`

  * These custom databases will be used to perform additional annatation processes using BLASTn

.. note::

   Users must must carefully read the documentation in order to better understand the details of the pipeline workflow customization.
