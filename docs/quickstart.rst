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

Prepare a samplesheet
---------------------

After downloading the genome, we must create a samplesheet for the input data as described in the :ref:`samplesheet manual page<samplesheet>`. A proper formated file for this data would look like that:

.. code-block:: yaml

  samplesheet:
    - id: ecoli
      assembly: ecoli_ref.fna
      resfinder: Escherichia coli

.. note::

  Download this file and save it as ``bacannot_samplesheet.yaml``.

Download databases
------------------

.. code-block:: bash

  # Download pipeline databases
  nextflow run fmalmeida/bacannot \
    --get_dbs \
    --output bacannot_dbs \
    -profile docker

.. note::

  You can choose between docker and singularity profiles. These are the possibilities right now. Conda shall come in near future.

Run the pipeline
----------------

In this step we will get a major overview of the main pipeline's steps. To run it, we will use the databases (``bacannot_dbs``) downloaded in the previous step.

.. code-block:: bash

  # Run the pipeline using the Escherichia coli resfinder database
  nextflow run fmalmeida/bacannot \
    --input bacannot_samplesheet.yaml \
    --output _ANNOTATION \
    --bacannot_db ./bacannot_dbs \
    --max_cpus 10

.. note::

  The resfinder species could also be selected via the command line with ``--resfinder_species``. Please, read more about it at :ref:`manual` and :ref:`samplesheet`.

Outputs
-------

A glimpse over the main outputs produced by bacannot is given at :ref:`outputs` section.

Testing more workflows
----------------------

Moreover, we have also made available a few example datasets in the pipeline so users can test all capabilities at once, from assembling raw reads to annotating genomes. To test it users must run:

.. code-block:: bash

  # Run the pipeline using the provided (bigger) test dataset
  nextflow run fmalmeida/bacannot --profile docker,test --bacannot_db ./bacannot_dbs --max_cpus 10

  # Or run the quick test
  nextflow run fmalmeida/bacannot --profile docker,quicktest --bacannot_db ./bacannot_dbs ---max_cpus 10

.. note::

  Unfortunately, due to file sizes, we could not provide fast5 files for users to check on the methylation step.