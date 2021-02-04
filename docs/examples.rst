.. _examples:

Usage examples
==============

Launching interactive graphical interface
"""""""""""""""""""""""""""""""""""""""""

Users can trigger a graphical and interactive pipeline configuration and execution by using `nf-core launch <https://nf-co.re/launch>`_ utility.

.. code-block:: bash

     # Install nf-core
     pip install nf-core

     # Lauch pipeline interactive configuration
     nf-core launch fmalmeida/bacannot

Single genome annotation
""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/bacannot --outdir TESTE --threads 3 --genome assembly.fasta \
      --bedtools_merge_distance -20 --skip_kofamscan

.. note::

  This command will perform a rapid annotation of ``assembly.fasta`` file using a minimum of 20 overlapping bases
  for gene merge and will not execute Kofamscan, nor methylation call with Nanopolish.

Multiple genome annotation
""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/bacannot --outdir TESTE --threads 3 --in_yaml samplesheet.yaml \
      --custom_db db1.fasta

.. warning::

  Samplesheet must be properly configured as in :ref:`samplesheet`.

.. note::

  The ``--custom_db`` parameter is used to add an annotation process with BLASTn using an user's custom db.

A little more complex example
"""""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/bacannot --outdir TESTE --threads 3 --genome assembly.fasta --bedtools_merge_distance -20 \
      --nanopolish_fastq_reads "fastq/input.fastq" --nanopolish_fast5_dir "fast5_pass_dir" --resfinder_species "Escherichia coli"

.. note::

  Differently, this command will run **all** the main analysis because the Resfinder and Nanopolish
  parameters have been set and no process have been told to skip (e.g. ``--skip_kofamscan``).

Annotating from raw reads
"""""""""""""""""""""""""

Users are able to annotate genomes directly from raw reads. When raw reads are used, Unicycler is used to create
shortreads-only and hybrid assemblies while Flye is used to create longreads-only assemblies the annotation process.


::

      nextflow run fmalmeida/bacannot --sreads_paired "sample1_{1,2}.fastq" --lreads "sample1_lreads.fastq" --lreads_type nanopore \
      --outdir TESTE --skip_kofamscan --threads 5 --nanopolish_fastq_reads "sample1_lreads.fastq" --nanopolish_fast5_dir "fast5_pass_dir"

.. note::

  This command will first perform a hybrid assembly with Unicycler and then annotate the assembled genome. Additionnally, since
  nanopolish parameters were given, it will call methylations with nanopolish.

Running with a configuration file
"""""""""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/bacannot -c bacannot.config
