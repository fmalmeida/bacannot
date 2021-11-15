.. _samplesheet:

Samplesheet (input files)
=========================

The samplesheet is a required YAML document that is used to describe the input samples and, if desired, its "sample-specific" configuration. The input samplesheet is given using the ``--input`` parameter.

.. tip::

  A samplesheet template can be downloaded with: nextflow run fmalmeida/bacannot --get_samplesheet

Samplesheet header
""""""""""""""""""

The first line of the file must be the header followed by an indentation:

.. code-block:: yaml

  samplesheet:
    - ...:

Sample identification
"""""""""""""""""""""

Each sample must be identified by the tag *id* in the YAML file, followed by the input input tags that shall
be used by the pipeline:

.. code-block:: yaml

  samplesheet:
    - id: sample_1
      ...:
      ...:
    - id: sample_2
      ...:
      ...:

Input tags (keys)
"""""""""""""""""

Input tags are are used to represent/set the inputs that shall be used for each input sample. By default, for resfinder species panel, if it is not set inside the samplesheet, the pipeline will use the configurations set via the "nextflow config file" or via the command line. Otherwise, if set inside the samplesheet, it will overwrite the pipeline’s configuration for that specific sample.

Please, the :ref:`manual reference page<manual>` the global/defaults configurations.

The available tags are:

.. list-table::
   :widths: 20 50
   :header-rows: 1

   * - Input tags
     - Description

   * - ``assembly``
     - Used to set path to genomic FASTA of an assembled bacterial genome

   * - ``illumina``
     - Used to set path to illumina raw reads (paired, unpaired or both)

   * - ``pacbio``
     - Used to set path to pacbio raw reads (mutually excludable with ``nanopore``)

   * - ``nanopore``
     - Used to set path to nanopore raw reads (mutually excludable with ``pacbio``)

   * - ``fast5``
     - Used to set path to nanopore raw FAST5 data (used in conjunction with ``nanopore`` for calling methylation with Nanopolish)

   * - ``resfinder``
     - Used to set resfinder species database for resistance annotation with resfinder (must be exactly as shown in `their web page <https://cge.cbs.dtu.dk/services/ResFinder/>`_). If your species is not available in Resfinder panels, you may use it with the "Other" panel.


.. note::

  Note for the illumina tag/key.

  * When using both paired and unpaired reads, the paired reads must be given first, in the order\: pair 1, pair 2, unpaired.
  * Otherwise, if using only paired reads, they must be given in the order\: pair 1, pair 2.
  * If using only unpaired reads, only one entry is expected. Check samples in the template to 1, 4 and 5 to understand it.
  * The illumina tag is the only one that **must** be set in indented newlines
      * two white spaces relative to the
      * one line per read as shown in the complete samplesheet example.

.. warning::

  All the other input tags **must** be set in the same line, right after the separator (":"), without quotations, white spaces or signs.

Complete samplesheet example
""""""""""""""""""""""""""""

.. code-block:: yaml

  samplesheet:
    - id: sample_1
      illumina:
        - sample_1/1.fastq
        - sample_1/2.fastq
      nanopore: sample_1/ont.fastq
    - id: sample_2
      assembly: sample_2/assembly.fasta
      nanopore: sample_2/ont.fastq
      fast5: sample_2/fast5_pass
      resfinder: Klebsiella              # this tells the pipeline a differente value for only this sample
    - id: sample_3
      nanopore: sample_3/ont.fastq
      fast5: sample_3/fast5_pass
    - id: sample_4
      pacbio: sample_4/pacbio.fastq
      illumina:
        - sample_4/merged_unpaired.fastq
    - id: sample_5
      illumina:
        - sample_5/1.fastq
        - sample_5/2.fastq
        - sample_5/merged.fastq
