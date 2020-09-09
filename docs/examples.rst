.. _examples:

CLI usage Examples
==================

Simple annotation example
"""""""""""""""""""""""""

::

      ./nextflow run main.nf --outdir TESTE --threads 3 --genome assembly.fasta \
      --bedtools_merge_distance -20 --not_run_kofamscan

.. note::

  This command will perform a rapid annotation of ``assembly.fasta`` file using a minimum of 20 overlapping bases
  for gene merge and will not execute KofamScan, nor methylation call with Nanopolish.

A little more complex example
"""""""""""""""""""""""""""""

::

      ./nextflow run main.nf --outdir TESTE --threads 3 --genome assembly.fasta --bedtools_merge_distance -20 \
      --nanopolish_fastq_reads "fastq/input.fastq" --nanopolish_fast5_dir "fast5_pass_dir"

.. note::

  Differently, this command will run **all** analysis because the Nanopolish parameters have
  been set and no process have been told to skip (e.g. ``--not_run_kofamscan``).


Running with a configuration file
"""""""""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/bacannot -c bacannot.config
