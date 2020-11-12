.. _examples:

CLI usage Examples
==================

Simple annotation example
"""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/bacannot --outdir TESTE --threads 3 --genome assembly.fasta \
      --bedtools_merge_distance -20 --not_run_kofamscan

.. note::

  This command will perform a rapid annotation of ``assembly.fasta`` file using a minimum of 20 overlapping bases
  for gene merge and will not execute Kofamscan, nor methylation call with Nanopolish.

More than one genome at once
""""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/bacannot --outdir TESTE --threads 3 --genome "inputs/*.fasta" --resfinder_species "Klebsiella"

.. warning::

  This option is incompatible with the methylation calling step. When running the methylation analysis users **MUST**
  run the pipeline with **only one** sample at a time, to ensure the right combination of files.

.. tip::

  The use of ``--resfinder_species`` parameter will activate the resfinder annotation process

A little more complex example
"""""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/bacannot --outdir TESTE --threads 3 --genome assembly.fasta --bedtools_merge_distance -20 \
      --nanopolish_fastq_reads "fastq/input.fastq" --nanopolish_fast5_dir "fast5_pass_dir"

.. note::

  Differently, this command will run **all** analysis because the Nanopolish parameters have
  been set and no process have been told to skip (e.g. ``--not_run_kofamscan``).

.. warning::

  When running the methylation analysis users **MUST** run the pipeline with **only one** sample at a time,
  to ensure the right combination of files.

Annotating from raw reads
"""""""""""""""""""""""""

Users are able to annotate genomes directly from raw reads. When raw reads are used, Unicycler is used to create
shortreads-only and hybrid assemblies while Flye is used to create longreads-only assemblies the annotation process.


::

      nextflow run fmalmeida/bacannot --sreads_paired "sample1_{1,2}.fastq" --lreads "sample1_lreads.fastq" --lreads_type nanopore \
      --outdir TESTE --not_run_kofamscan --threads 5 --nanopolish_fastq_reads "sample1_lreads.fastq" --nanopolish_fast5_dir "fast5_pass_dir"

.. warning::

  When combining different sequencing library types users **MUST** run the pipeline with **only one** sample at a time,
  to ensure the right combination of files.

.. warning::

  When running the methylation analysis users **MUST** run the pipeline with **only one** sample at a time,
  to ensure the right combination of files.

Running with an interactive graphical interface
"""""""""""""""""""""""""""""""""""""""""""""""

::

     nf-core launch fmalmeida/bacannot


Running with a configuration file
"""""""""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/bacannot -c bacannot.config
