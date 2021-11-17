.. _installation:

Installation
============

Dependencies
------------

This pipeline requires only `Docker <https://www.docker.com/>`_ (and its Docker images) and
`Nextflow <https://www.nextflow.io/docs/latest/index.html>`_ to run.

1. Installing Docker

   * Read more in their `manual <https://docs.docker.com/>`_

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

      docker pull fmalmeida/bacannot:main_tools ;  # this is the core of the main image
      docker pull fmalmeida/bacannot:v3.0       ;
      docker pull fmalmeida/bacannot:kofamscan  ;
      docker pull fmalmeida/bacannot:antismash  ;
      docker pull fmalmeida/bacannot:jbrowse    ;
      docker pull fmalmeida/bacannot:v3.0_renv  ;

.. tip::

   If the download of ``fmalmeida/bacannot:v3.0`` image keeps hanging due to its size, download the ``fmalmeida/bacannot:main_tools`` first. It is the core of the versioned tag and it will help on the download by creating some cache. Also, remember to always keep your Docker images up to date (Docker pull will always download the latest)

6. (Optional) Install nf-core utility

   .. code-block:: bash

      pip install nf-core>=1.10

.. note::

  The pipeline requires a UNIX system, therefore, Windows users may successfully use this pipeline via the `Linux subsystem for windows <https://docs.microsoft.com/pt-br/windows/wsl/install-win10>`_.

  Nextflow team has made available a `nice tutorial <https://www.nextflow.io/blog.html>`_ about this issue.
