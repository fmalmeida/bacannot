/*
 * Define help message
 */

def helpMessage() {
  log.info """
  Usage:
  nextflow run fmalmeida/bacannot [--help] [ -c nextflow.config ] [OPTIONS] [-with-report] [-with-trace] [-with-timeline]
   
  The command line help message is not supported anymore because it has become too big and, therefore, very clumsy and confusing to read.

  Please, use the nextflow.config configuration file that has comments and help messages inside it or refer to the online manual:
    
            https://bacannot.readthedocs.io/en/latest/manual.html
    
  Comments:
  This pipeline contains a massive amount of configuration variables and its usage as CLI parameters could cause the command line to be huge.
  
  Therefore, it is recommended to use the nextflow.config file in order to make parameterization easier and more readable.

  To create a template configuration or samplesheet file, use:
  nextflow run fmalmeida/bacannot [--get_config] [--get_samplesheet]
  """.stripIndent()
}
