.. _samplesheet:

Samplesheet configuration (for multi-genome analysis)
=====================================================

The samplesheet is a YAML document that is used to describe the input samples. It is required when the user
wants to annotate more than one genome at once saving them at the same output directory. This execution is
triggered by the ``--in_yaml`` parameter and it is incompatible with all the parameters used for single
genome analysis (shown below):

The use of a samplesheet is incompatible with:

* --genome
* --sreads_paired
* --sreads_single
* --lreads
* --lreads_type
* --resfinder_species
* --nanopolish_fast5
* --nanopolish_fastq

Samplesheet example
"""""""""""""""""""

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
