.. _custom-db:

Custom database configuration
=============================

It is also possible that users use custom databases for the annotation of genomes. Currently, the pipeline only executes BLASTn alignments using the genome as query.
Therefore, the given databases must be files of **target gene sequences in nucleotide** FASTA.

The custom annotation is triggered with the ``--custom_db`` parameter. The pipeline accepts more than one custom database at once, separated by commas, e.g.
``--custom_db db1.fasta,db2.fasta``.

Although simple, the custom database must follow some rules about sequence header format in order to make it possible the summarization of alignments and renderization
of custom reports in HTML format, that shall be available under the ``report_files`` directory.

Sequence header format
""""""""""""""""""""""
