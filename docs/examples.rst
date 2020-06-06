.. _examples:

CLI usage Examples
==================

Simple annotation example
"""""""""""""""""""""""""

::

      ./nextflow run main.nf --outDir TESTE --threads 3 --genome assembly.fasta \
      --bedtools_merge_distance -20 --not_run_kofamscan

.. note::

  This command will perform a rapid annotation of ``assembly.fasta`` file using a minimum of 20 overlapping bases
  for gene merge and will not execute KofamScan, nor methylation call with Nanopolish, nor pangenome analysis with
  Roary since their required arguments are not defined.

A little more complex example
"""""""""""""""""""""""""""""

::

      ./nextflow run main.nf --outDir TESTE --threads 3 --genome assembly.fasta --bedtools_merge_distance -20 \
      --roary_reference_genomes "references/*.fasta" --nanopolish_fastq_reads "fastq/input.fastq" \
      --nanopolish_fast5_dir "fast5_pass_dir" --diamond_minimum_alignment_length 500

.. note::

  Differently, this command will run **all** analysis because all optional arguments of Roary and Nanopolish have
  been set and no process have been told to skip (e.g. ``--not_run_kofamscan``). Additionally, we have set a new
  minimum diamond alignment length to report hits (500).


Running with a configuration file
"""""""""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/bacannot -c bacannot.config
