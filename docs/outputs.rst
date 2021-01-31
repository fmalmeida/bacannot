.. _outputs:

Outputs overview
================

Following the same results produced in the :ref:`quickstart` section, the outputs are presented again here to provide a specific page for them. The quickstart
command wrote the results under the ``_ANNOTATION`` directory.

.. note::

  Please take note that the pipeline uses the directory set with the ``--outdir`` parameter as a storage place in which it will create a folder named as the
  ``--prefix`` parameter. This ``{prefix}`` folder will contain all the results. Therefore the the same ``--outdir`` can be used for different runs|genomes
  as each one of them will have a different sub-folder. This is useful and required for the genomic comparative pipeline (that is under construction) that will
  use this folder as input, and enable the user to rapidly compare the results between the samples under the same ``--outdir`` folder.

Directory tree
--------------

After a successful execution, you will have something like this:

.. code-block:: bash

    # Directory tree from the running dir
    .
    ├── _ANNOTATION
    │   └── ecoli
    │       ├── ICEs                                            # Results from ICEberg database annotation
    │       ├── KOfamscan                                       # Results from annotation with KEGG database
    │       ├── MLST                                            # MLST results with mlst pipeline
    │       ├── annotation                                      # Prokka annotation files
    │       ├── assembly                                        # Assembly files (when raw reads are given)
    │       ├── gbk                                             # Gbk file produced from the resulting GFF
    │       ├── genomic_islands                                 # Genomic Islands predicted with IslandPath-DIMOB
    │       ├── gffs                                            # A copy of the main GFF files produced during the annotation
    │       ├── jbrowse                                         # The files that set up the JBrowse genome browser
    │       ├── plasmids                                        # Plasmid annotation results from Platon and Plasmidfinder
    │       ├── prophages                                       # Prophage annotation results from PhiSpy, Phigaro and PHAST
    │       ├── rRNA                                            # barrnap annotation results
    │       ├── report_files                                    # Annotation reports in HTML format
    │       ├── resistance                                      # AMR annotation results from ARGminer, AMRFinderPlus, RGI and Resfinder
    │       ├── sqldb                                           # The sqlDB of the annotation used by the shiny server for rapid parsing
    │       ├── tools_versioning                                # Versions of tools and databases used (whenever available)
    │       ├── virulence                                       # Virulence genes annotation results from Victors and VFDB databases
    │       └── run_server.sh                                   # The shiny parser runner that enables a rapid and simple exploration of the results (see below)
    └── oxford.fasta

Bacannot shiny parser
---------------------

The bacannot shiny server is basically a wrapper of the main outputs of the pipeline that is packed up in a docker image called ``fmalmeida/bacannot:server``.
This server is triggered by going under the results folder, in our quickstart case, the ``_ANNOTATION/ecoli`` folder, and executing the command:

.. code-block:: bash

  # Trigger the server
  ./run_server.sh -s

  # This will open the pipeline in localhost:3838
  # log message:
  The server has started in: http://localhost:3838/
  When finished, run the command:
	       docker rm -f {docker container id}

  # To stop the server you just need to execute
  docker rm -f {docker container id}

Server homepage
^^^^^^^^^^^^^^^

In the first page it has indexed as url links the main HTML reports and the JBrowse genome browser.

.. image:: images/bacannot_server_home.png
  :width: 800
  :align: center

Server sqlDB parser
^^^^^^^^^^^^^^^^^^^

In the second page, the sqlDB is used to provide a rapid and simple way to query and filter the genome annotation.

.. note::

  The sqlDB parser contains a set of features that enables that the users filter the annotation following their desires. It is possible
  to filter based on the ``contigs``, ``sources``, ``start``, ``end``, ``strand`` and more.

  Additionally, the parser accepts as input a file of patterns to filter the annotation based on the values available in the attributes
  column of the GFF (9th column). Any value available in this column can be used as filters, the only requirement is to write each pattern
  in one line, exactly as it is found in the annotation result. For example, it can be used to select only a few genes based on their IDs.


.. image:: images/bacannot_server_sqldb.png
  :width: 800
  :align: center


Server BLAST app
^^^^^^^^^^^^^^^^

In the last page, the server provides a simple way to BLAST the genome with new gene queries and to automatically identify intersections
between the blast results and the the main annotation.

.. image:: images/bacannot_server_blast.png
  :width: 800
  :align: center
