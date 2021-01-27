.. _quickstart:

Quickstart
==========

For a rapid and complete quickstart we will use as input the nanopore raw reads provided in the `Canu quickstart section <https://canu.readthedocs.io/en/latest/quick-start.html#assembling-pacbio-clr-or-nanopore-data>`_.

Download the data
"""""""""""""""""

.. code-block:: bash

  # Download and save as oxford.fasta
  curl -L -o oxford.fasta http://nanopore.s3.climb.ac.uk/MAP006-PCR-1_2D_pass.fasta

Run the pipeline
""""""""""""""""

For examplification purposes and to get a major overview we will execute the pipeline's major processes:

.. code-block:: bash

  # Run the pipeline using the Escherichia coli resfinder database
  nextflow run fmalmeida/bacannot --prefix example \
  --lreads oxford.fasta \
  --lreads_type nanopore \
  --outdir example_out \
  --threads 4 \
  --resfinder_species "Escherichia coli"

Outputs
"""""""

This command shall write the results under the ``example_out`` directory.
