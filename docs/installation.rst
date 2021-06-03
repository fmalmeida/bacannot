.. _installation:

Installation
============

Dependencies
------------

This pipeline requires only `Docker <https://www.docker.com/>`_ (and its Docker images) and
`Nextflow <https://www.nextflow.io/docs/latest/index.html>`_ to run.

1. Installing Docker

  + Read more in their `manual <https://docs.docker.com/>`_

2. Installing Nextflow

    .. code-block:: bash

      curl -s https://get.nextflow.io | bash

3. Download the pipeline

    .. code-block:: bash

      nextflow pull fmalmeida/bacannot

4. Test your installation

    .. code-block:: bash

      nextflow run fmalmeida/bacannot --help

5. Download required Docker images

    .. code-block:: bash

      docker pull fmalmeida/bacannot:v2.3 ;
      docker pull fmalmeida/bacannot:kofamscan ;
      docker pull fmalmeida/bacannot:jbrowse ;
      docker pull fmalmeida/bacannot:v2.3_renv

6. (Optional) Install nf-core utility

    .. code-block:: bash

      pip install nf-core>=1.10

7. (Optional) Docker image for using raw reads as input

    .. code-block:: bash

      docker pull fmalmeida/mpgap

.. note::

  Now, everything is set up and ready to run. Remember to always keep your Docker images up to date (Docker pull will always download the latest).

  The pipeline requires a UNIX system, therefore, Windows users may successfully use this pipeline via the `Linux subsystem for windows <https://docs.microsoft.com/pt-br/windows/wsl/install-win10>`_.
