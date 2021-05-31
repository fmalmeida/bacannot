.. _samplesheet:

Samplesheet configuration (for multi-genome analysis)
=====================================================

The samplesheet is a YAML document that is used to describe the input samples. It is required when the user
wants to annotate more than one genome at once saving them at the same output directory. This execution is
triggered by the ``--in_yaml`` parameter and it is incompatible with all the parameters used for single
genome analysis (shown below):

The use of a samplesheet is **incompatible** with:

+ ``--genome``
+ ``--sreads_paired``
+ ``--sreads_single``
+ ``--lreads``
+ ``--lreads_type``
+ ``--resfinder_species``
+ ``--nanopolish_fast5``
+ ``--nanopolish_fastq``

Samplesheet header
""""""""""""""""""

The first line of the file must be the header followed by an indentation:

.. code-block:: yaml

  samplesheet:
    - ...:

Sample identification
"""""""""""""""""""""

Each sample must be identified by the tag *id* in the YAML file, followed by the input file tags that shall
be used by the pipeline:

.. code-block:: yaml

  samplesheet:
    - id: sample_1
      ...:
      ...:
    - id: sample_2
      ...:
      ...:

File tags
"""""""""

File tags are the tags that are used to represent/set the input files that shall be used for each sample that
will be analysed. The available file tags are:

.. list-table::
   :widths: 20 50
   :header-rows: 1

   * - File tags
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
     - Used to set resfinder species database for resistance annotation with resfinder (must be exactly as shown in their manual/web tool)


.. note::

  The illumina tag is the only one that **must** be set in indented newlines (one line per read) as shown in the complete samplesheet example. The order
  of the reads in these newlines must be Pair1; Pair2; Unpaired (Whenever they are used) -- Check samples 1, 4 and 5 to understand.

  All the other file tags **must** be set in the same line, right after the separator (":"), without quotations.

Complete samplesheet example
""""""""""""""""""""""""""""

.. code-block:: yaml

  samplesheet:
    - id: sample_1
      illumina:
        - sample_1/1.fastq
        - sample_1/2.fastq
      nanopore: sample_1/ont.fastq
      resfinder: Escherichia coli
    - id: sample_2
      assembly: sample_2/assembly.fasta
      nanopore: sample_2/ont.fastq
      fast5: sample_2/fast5_pass
      resfinder: Klebsiella
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
