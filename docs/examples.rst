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

.. code-block:: bash

      ./nextflow run fmalmeida/bacannot \
        --outdir TESTE \
        --threads 3 \
        --genome assembly.fasta \
        --bedtools_merge_distance -20 \
        --skip_kofamscan

.. note::

  This command will perform a rapid annotation of ``assembly.fasta`` file using a minimum of 20 overlapping bases
  for gene merge and will not execute Kofamscan, nor methylation call with Nanopolish.

Multiple genome annotation
""""""""""""""""""""""""""

.. code-block:: bash

      ./nextflow run fmalmeida/bacannot \
        --outdir TESTE \
        --threads 3 \
        --in_yaml samplesheet.yaml \
        --custom_db db1.fasta

.. warning::

  Samplesheet must be properly configured as in :ref:`samplesheet`.

.. note::

  The ``--custom_db`` parameter is used to add an annotation process with BLASTn using an user's custom db.

A little more complex example
"""""""""""""""""""""""""""""

.. code-block:: bash

      ./nextflow run fmalmeida/bacannot \
        --outdir TESTE \
        --threads 3 \
        --genome assembly.fasta \
        --bedtools_merge_distance -20 \
        --nanopolish_fastq "fastq/input.fastq" \
        --nanopolish_fast5 "fast5_pass_dir" \
        --resfinder_species "Escherichia coli"

.. note::

  Differently, this command will run **all** the main analysis because the Resfinder and Nanopolish
  parameters have been set and no process have been told to skip (e.g. ``--skip_kofamscan``).

Annotating from raw reads
"""""""""""""""""""""""""

Users are able to annotate genomes directly from raw reads. When raw reads are used, Unicycler is used to create
shortreads-only and hybrid assemblies while Flye is used to create longreads-only assemblies the annotation process.


.. code-block:: bash

      nextflow run fmalmeida/bacannot \
        --sreads_paired "sample1_{1,2}.fastq" \
        --lreads "sample1_lreads.fastq" \
        --lreads_type nanopore \
        --outdir TESTE \
        --skip_kofamscan \
        --threads 5 \
        --nanopolish_fastq "sample1_lreads.fastq" \
        --nanopolish_fast5 "fast5_pass_dir"

.. note::

  This command will first perform a hybrid assembly with Unicycler and then annotate the assembled genome. Additionnally, since
  nanopolish parameters were given, it will call methylations with nanopolish.

.. note::
  
  Remember to always write input paths inside double quotes.

.. note::
  
  When using paired end reads it is required that input reads are set with the "{1,2}"" pattern. For example: "SRR6307304_{1,2}.fastq". This will properly load reads "SRR6307304_1.fastq" and "SRR6307304_2.fastq"

.. warning::
  
  When running hybrid assemblies or mixing short read types it is advised to **avoid not required REGEX** and write the full file path, using only the required REGEX for paired end reads when applicable. So that the pipeline does not load any different read that also matches the REGEX and avoid confusions with the inputs.

Running with a configuration file
"""""""""""""""""""""""""""""""""

.. code-block:: bash

      ./nextflow run fmalmeida/bacannot -c bacannot.config
