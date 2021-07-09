.. _quickstart:

Quickstart
==========

For a rapid and simple quickstart we will use as input the *Escherichia coli* reference genome.

Download the data
-----------------

.. code-block:: bash

  # Download the ecoli ref genome
  wget -O ecoli_ref.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz
  gzip -d ecoli_ref.fna.gz

Run the pipeline
----------------

For examplification purposes and to get a major overview we will execute the pipeline's major processes:

.. code-block:: bash

  # Run the pipeline using the Escherichia coli resfinder database
  nextflow run fmalmeida/bacannot \
    --prefix ecoli \
    --genome ecoli_ref.fna \
    --outdir _ANNOTATION \
    --threads 4 \
    --resfinder_species "Escherichia coli"

Outputs
-------

A glimpse over the main outputs produced by bacannot is given at :ref:`outputs` section.
