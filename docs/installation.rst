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

5. Download pipeline databases

   .. code-block:: bash

      nextflow run fmalmeida/bacannot -profile docker --get_dbs

.. tip::

   This will download the required databases and save them in the current working directory. To save it somewhere else, use ``--output``.

6. Download required Docker images

   .. code-block:: bash

      docker pull fmalmeida/bacannot:v3.1_misc    ;
      docker pull fmalmeida/bacannot:v3.1_perlenv ;
      docker pull fmalmeida/bacannot:v3.1_pyenv   ;
      docker pull fmalmeida/bacannot:v3.1_py36env ;
      docker pull fmalmeida/bacannot:v3.1_renv    ;
      docker pull fmalmeida/bacannot:jbrowse      ;

.. tip::

   Required images can be downloaded on the fly, not requiring to be previously available. Just be sure to select the correct desired profile (``docker/singularity``).

.. note::

  The pipeline requires a UNIX system, therefore, Windows users may successfully use this pipeline via the `Linux subsystem for windows <https://docs.microsoft.com/pt-br/windows/wsl/install-win10>`_.

  Nextflow team has made available a `nice tutorial <https://www.nextflow.io/blog.html>`_ about this issue.
