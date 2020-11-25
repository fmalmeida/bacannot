def paramsCheck() {
  // Input parameters
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

  // Checking the use of lreads
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

  // Prokka parameters
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

  // Methylation parameters
  if ((params.nanopolish_fast5_dir && !params.nanopolish_fastq_reads) || (!params.nanopolish_fast5_dir && params.nanopolish_fastq_reads)) {
    println """
    ERROR!

    A minor error has occurred
      ==> User have forget to set both --nanopolish_fast5_dir and --nanopolish_fastq_reads.

    These parameters must be used together. They are the necessary files to call methylations from ONT data with Nanopolish.

    Cheers.
    """.stripIndent()

    exit 1
  }

  params.help = false
   // Show help emssage
   if (params.help){
     helpMessage()
     exit 0
  }

  // CLI examples
  params.examples = false
   // Show help emssage
   if (params.examples){
     exampleMessage()
     exit 0
  }

  /*
   * Does the user wants to download the configuration file?
   */

  params.get_config = false
  if (params.get_config) {
    new File("bacannot.config").write(new URL ("https://github.com/fmalmeida/bacannot/raw/master/nextflow.config").getText())
    println ""
    println "bacannot.config file saved in working directory"
    println "After configuration, run:"
    println "nextflow run fmalmeida/bacannot -c ./bacannot.config"
    println "Nice code!\n"
    exit 0
  }
}
