.. _installation:

Installation
============

Dependencies
------------

This pipeline requires only `Docker <https://www.docker.com/>`_ (and its Docker images) and
`Nextflow <https://www.nextflow.io/docs/latest/index.html>`_ to run.

1. Installing Docker
  * Read more in their `manual <https://docs.docker.com/>`_
  * Or give this `in-house script <https://github.com/fmalmeida/bioinfo/blob/master/dockerfiles/docker_install.sh>`_ a try.
2. Installing Nextflow

    ``curl -s https://get.nextflow.io | bash``

3. Download the pipeline

    ``./nextflow pull fmalmeida/bacannot``

4. Test your installation

    ``./nextflow run fmalmeida/bacannot --help``

5. Download required Docker images

    ``docker pull fmalmeida/compgen:BACANNOT``

    ``docker pull fmalmeida/compgen:KOFAMSCAN``

    ``docker pull fmalmeida/compgen:JBROWSE``

    ``docker pull fmalmeida/compgen:RENV``

Now, everything is set up and ready to run. Remember to always keep your Docker images up to date (Docker pull you always download the latest).
