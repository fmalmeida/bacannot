.. _input:

Input files
===========

Users can perform the annotation analysis using either raw reads or assembled genomes as input. When raw reads are used, Unicycler is used to create
shortreads-only and hybrid assemblies while Flye is used to create longreads-only assemblies the annotation process.

* path to genome fasta file **OR** to raw reads.
* path to a directory of FAST5 and path to ONT fastq to be used for methylation calling (optional)
* In order to get the best results from this pipeline, users are advised to analyse one sample at a time.

.. note::

  Users can analyse more than one genome at once by configuring a samplesheet and
  giving it as input with the ``--in_yaml`` parameter. See :ref:`samplesheet`.

.. note::

   Users must **never** use hard or symbolic links. This will make nextflow fail.
   When setting the parameters, please **always** give full path to a hard file,
   not to a link. This will prevent file access fail.
