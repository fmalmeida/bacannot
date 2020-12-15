def paramsCheck() {

  /*

      User tried to use both raw reads and a assembled genome as input

  */
  if (params.genome && (params.sreads_paired || params.sreads_single || params.lreads)) {
    println """
    ERROR!

    A minor error has occurred
      ==> User used raw reads and assembled genomes as input.

    You cannot use both types of inputs together. Used either assembled genomes OR raw reads.

    Cheers.
    """.stripIndent()

    exit 1
  }

  /*

      User tried to use both the options for single and multiple genome analysis

  */
  if (params.in_yaml && (params.sreads_paired || params.sreads_single || params.lreads
                         || params.genome || params.lreads_type || params.resfinder_species
                         || params.nanopolish_fast5 || params.nanopolish_fastq)) {
    println """
    ERROR!

    A major error has occurred
      ==> User have set parameters for single and multiple genome analysis.

    This pipeline works either annotating a single genome or multiple genome at once. However, the analysis of multiple
    genomes is set via the YAML file and it is incompatible with the parameters from single genome analysis.

    Therefore, when using the --in_yaml parameter you cannot use any of the following parameters as they are specific for single genome analysis:

      * --genome
      * --sreads_paired
      * --sreads_single
      * --lreads
      * --lreads_type
      * --resfinder_species
      * --nanopolish_fast5
      * --nanopolish_fastq

    Cheers.
    """.stripIndent()

    exit 1
  }

  /*

      User has given the --lreads parameter but forgot --lreads_type

  */
  if (params.lreads && !params.lreads_type) {
    println """
    ERROR!

    A minor error has occurred
      ==> User used --lreads but forgot --lreads_type.

    When giving longreads as input, you must tell the pipeline from wich tech it comes from: 'nanopore' or 'pacbio'

    Cheers.
    """.stripIndent()

    exit 1
  }

  /*

      Checking the prokka parameters

  */
  if (params.prokka_kingdom && !params.prokka_genetic_code) {
    println """
    ERROR!

    A minor error has occurred
      ==> User have set --prokka_kingdom but forget --prokka_genetic_code.

    These parameters must be used together. If you change prokka defaults kingdom parameter you must set the genetic code to be used for translation.

    If in doubt with these parameters let it blank, or get more information in Prokka's documentation.

    Cheers.
    """.stripIndent()

    exit 1
  }

  /*

      Checking nanopolish parameters

  */
  if ((params.nanopolish_fast5 && !params.nanopolish_fastq) || (!params.nanopolish_fast5 && params.nanopolish_fastq)) {
    println """
    ERROR!

    A minor error has occurred
      ==> User forgot to set both --nanopolish_fast5 and --nanopolish_fastq.

    These parameters must be used together. They are the necessary files to call methylations from ONT data with Nanopolish.

    Cheers.
    """.stripIndent()

    exit 1
  }

  /*

      Checking resfinder parameters

  */
  if (params.resfinder_species == "other" || params.resfinder_species == "Other" || params.resfinder_species == "OTHER") {
    println """
    ERROR!

    A minor error has occurred
      ==> User has set the resfinder panel to "Other"

    This is impossible, since the pipeline tries to annotation point finder mutation and these are incompatible with the "Other" panel

    Cheers.
    """.stripIndent()

    exit 1
  }

}
