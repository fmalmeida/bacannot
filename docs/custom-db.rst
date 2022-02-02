.. _custom-db:

Custom database configuration
=============================

It is also possible that users use custom databases for the annotation of genomes, either by providing custom nucleotide FASTAs, properly formatted with **target gene sequences in nucleotide** to be search against the genome with BLASTn, or by providing a file containing a list of NCBI Protein database accession which will retrieve these sequences and annotated the genome against them using BLASTp.

NCBI Protein database
---------------------

To perform an additional annotation of the genome using proteins from NCBI Protein database, one must use the parameter ``--ncbi_proteins`` and provide a file containing a list of protein IDs such as this one:

.. code-block:: bash

    WP_118891437.1
    VTX70803.1
    VTX49335.1
    WP_005693332.1
    WP_138172127.1

Custom nucleotide database (in FASTA)
-------------------------------------

This custom annotation is triggered with the ``--custom_db`` parameter. The pipeline accepts more than one custom database at once, separated by commas, e.g.
``--custom_db db1.fasta,db2.fasta``.

Although simple, the custom database must follow some rules about sequence header format in order to make it possible the summarization of alignments and renderization
of custom reports in HTML format, that shall be available under the ``report_files`` directory.

Sequence header format
""""""""""""""""""""""

Sequence headers must follow a 5-field rule separated by "~~~" and spaces. The first 4 fields must be separated by "~~~" and the last one by one space, following the
example shown below:

.. warning::

  Except for the Description field, the first four (DB name, gene name, gene reference and gene product) must be written **whithout** whitespaces.

.. code-block:: bash

    # Sequence header
    >Database_name~~~gene/product_name/alias~~~acc.number_reference_or_identification~~~gene_product  Description

    # An example with a VFDB sequence
    >VFDB~~~(plc)~~~VFG037176(gb|YP_001844723)~~~[Phospholipase_C_(VF0470)]  VFG037176(gb|YP_001844723) (plc) phospholipase C [Phospholipase C (VF0470)]

.. note::

  It is very important to follow this header format in order to make it possible and easier to render summaries and reports of the BLASTn result, such as below:

Summary example for custom database annotations
-----------------------------------------------

The custom annotations, either with ``--custom_db`` or with ``--ncbi_protein`` will produce useful summaries and reports (:ref:`outputs`) such as in this example:

.. code-block:: bash

    SEQUENCE	START	END	STRAND	GENE	COVERAGE	GAPS	%COVERAGE	%IDENTITY	DATABASE	VFDB_ID	PRODUCT	DESCRIPTION
    BNHFPCLP_00473	1	528	+	(rcsB)	1-527/651	2/3	80.8	99.43	VFDB	VFG049018(gb|YP_002920501.1)	[RcsAB_(VF0571)]	 VFG049018(gb|YP_002920501.1) (rcsB) transcriptional regulator RcsB [RcsAB (VF0571)]
    BNHFPCLP_00654	1	629	+	(cpsACP)	1-629/630	0/0	99.84	95.55	VFDB	VFG048985(gb|YP_002920368.1)	[Capsule_(VF0560)]	 VFG048985(gb|YP_002920368.1) (cpsACP) phosphatase PAP2 family protein [Capsule (VF0560)]
    BNHFPCLP_00676	1	1167	+	(ugd)	1-1167/1167	0/0	100.0	96.92	VFDB	VFG048797(gb|YP_002920350.1)	[Capsule_(VF0560)]	 VFG048797(gb|YP_002920350.1) (ugd) UDP-glucose 6-dehydrogenase [Capsule (VF0560)]
    BNHFPCLP_00680	1	741	+	(wzt)	1-741/741	0/0	100.0	99.46	VFDB	VFG049084(gb|YP_002920347.1)	[LPS_(VF0561)]	 VFG049084(gb|YP_002920347.1) (wzt) lipopolysaccharide O-antigen ABC transport system ATP-binding component [LPS (VF0561)]
    BNHFPCLP_00685	1	894	+	(wbbN)	1-894/894	0/0	100.0	98.77	VFDB	VFG049051(gb|YP_002920344.1)	[LPS_(VF0561)]	 VFG049051(gb|YP_002920344.1) (wbbN) glycosyltransferase [LPS (VF0561)]
    BNHFPCLP_00686	1	1131	+	(wbbO)	1-1131/1131	0/0	100.0	99.47	VFDB	VFG049040(gb|YP_002920343.1)	[LPS_(VF0561)]	 VFG049040(gb|YP_002920343.1) (wbbO) glycosyltransferase family 1 protein [LPS (VF0561)]
    BNHFPCLP_00770	1	624	+	(rcsA)	1-624/624	0/0	100.0	100.0	VFDB	VFG049007(gb|YP_002920216.1)	[RcsAB_(VF0571)]	 VFG049007(gb|YP_002920216.1) (rcsA) transcriptional activator for ctr capsule biosynthesis [RcsAB (VF0571)]
    BNHFPCLP_01728	1	936	+	(iroE)	1-936/936	0/0	100.0	99.04	VFDB	VFG044322(gb|YP_002919453)	[Sal_(VF0563)]	 VFG044322(gb|YP_002919453) (iroE) siderophore esterase IroE [Sal (VF0563)]
    BNHFPCLP_02092	1	543	+	(sciN/tssJ)	1-543/543	0/0	100.0	99.45	VFDB	VFG048784(gb|YP_005226619.1)	[T6SS_(VF0569)]	 VFG048784(gb|YP_005226619.1) (sciN/tssJ) type VI secretion system lipoprotein TssJ [T6SS (VF0569)]
