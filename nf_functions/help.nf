/*
 * Define help message
 */

def helpMessage() {
   log.info """
   Usage:
   nextflow run fmalmeida/bacannot [--help] [ -c nextflow.config ] [OPTIONS] [-with-report] [-with-trace] [-with-timeline]
   Comments:

   This pipeline contains a massive amount of configuration variables and its usage as CLI parameters would cause the command
   to be huge. Therefore, it is extremely recommended to use the nextflow.config configuration file in order to make
   parameterization easier and more readable.

   Creating a configuration file:

   nextflow run fmalmeida/bacannot [--get_config]

   Show command line examples:

   nextflow run fmalmeida/bacannot --examples

   Execution Reports:

   nextflow run fmalmeida/bacannot [ -c nextflow.config ] -with-report
   nextflow run fmalmeida/bacannot [ -c nextflow.config ] -with-trace
   nextflow run fmalmeida/bacannot [ -c nextflow.config ] -with-timeline

   OBS: These reports can also be enabled through the configuration file.
   OPTIONS:

      # Input configuration -- Analysis of a single genome
      # Users can give either a genome in FASTA file or raw reads in FASTQ
      # Please do not use glob. patterns ('*') with these parameters

    --prefix <string>                              Prefix for writing genome assembly and annotatin resulting files.
                                                   Preferentially the sample name. [Default: out]

    --genome <string>                              Set path to the genome in FASTA file.

    --sreads_paired <string>                       Illumina paired end reads in fastq for assembly before annotation.

    --sreads_single <string>                       Illumina unpaired reads in fastq for assembly before annotation.

    --lreads <string>                              Path to longreads in fastq assembly before annotation (ONT or Pacbio).

    --lreads_type <string>                         Tells the technology of the input longreads: [ nanopore or pacbio ].

      # Input configuration -- Analysis of multiple genomes
      # Users can give either a genome in FASTA file or raw reads in FASTQ
      # The analysis of multiple genomes at once is configured via a YAML file
      # Check the example YAML at: https://github.com/fmalmeida/bacannot/blob/master/example_samplesheet.yaml
      #
      # Also documented at: https://bacannot.readthedocs.io/en/latest/samplesheet.html

    --in_yaml <string>                             Set path to the samplesheet in YAML format to analyse more than one
                                                   genome at once.

      # Annotation configuration -- Used for either for the single
      # genome analysis workflow and the multiple genome analysis
      # Read it and configure it properly

      # General Parameters

    --outdir <string>                              Output directory name

    --threads <int>                                Number of threads to use

    --bedtools_merge_distance                      By default, this process is not executed. For execution
                                                   one needs to provide a value.Minimum number of overlapping
                                                   bases for gene merge using bedtools merge. Negative values,
                                                   such as -20, means the number of required overlapping bases
                                                   for merging. Positive values, such as 5, means the maximum
                                                   distance accepted between features for merging.


      # Prokka complementary parameters

    --prokka_kingdom <string>                      Prokka annotation mode. Possibilities (default 'Bacteria'):
                                                   Archaea|Bacteria|Mitochondria|Viruses

    --prokka_genetic_code <int>                    Genetic Translation code. Must be set if kingdom is not
                                                   default (in blank).

    --prokka_use_rnammer                           Tells prokka whether to use rnammer instead of barrnap.


      # Blast alignment parameters

    --blast_virulence_minid                        Min. identity % for virulence annotation. Default 90.

    --blast_virulence_mincov                       Min. gene/subject coverage for virulence annotation. Default 80.

    --blast_resistance_minid                       Min. identity % for resistance annotation. Default 90.

    --blast_resistance_mincov                      Min. gene/subject coverage for resistance annotation. Default 80.

    --blast_MGEs_minid                             Min. identity % for ICEs and prophage annotation. Default 65.

    --blast_MGEs_mincov                            Min. gene/subject coverage for ICEs and prophage annotation. Default 65.

    --plasmids_minid                               Min. identity % for plasmid detection. Default 90.

    --plasmids_mincov                              Min. coverage for plasmid detection. Default 60.

    --blast_custom_minid                           Min. identity % for the annotation using user's custom database. Default 0.

    --blast_custom_mincov                          Min. gene/subject coverage % for the annotation using user's custom database. Default 0.

      # User's custom database for annotation
      # Must be in gene nucleotide FASTA
      #
      # Well documented at: https://bacannot.readthedocs.io/en/latest/custom-db.html

    --custom_db                                    Path to the nucleotide FASTA file containing the user's custom database for annotation.
                                                   Multiple FASTAs can be provided separated by comma. E.g. db1.fasta,db2.fasta,...


      # Configure resfinder optional parameter
      # Only used with analysing a single genome
      # When analysing multiple genomes it must be set in the YAML file.
      # Check the example YAML at: https://github.com/fmalmeida/bacannot/blob/master/example_samplesheet.yaml
      #
      # Also documented at: https://bacannot.readthedocs.io/en/latest/samplesheet.html

    --resfinder_species                            It sets the species to be used for Resfinder annotation. If blank,
                                                   it will not be executed. Must be identical (without the *) as written
                                                   in their webservice https://cge.cbs.dtu.dk/services/ResFinder/.
                                                   E.g. 'Escherichia coli'; 'Klebsiella' ...

      # Configure (on/off) optional processes

    --skip_virulence_search                        Tells whether you do not want to execute virulence annotation

    --skip_resistance_search                       Tells whether you do not want to execute resistance annotation

    --skip_iceberg_search                          Tells whether you do not want to execute ICE annotation

    --skip_prophage_search                         Tells whether you do not want to execute prophage annotation

    --skip_plasmid_search                          Tells whether you do not want to execute plasmid detection

    --skip_kofamscan                               Tells whether you do not want to execute KO annotation with kofamscan


      # Configure optional Methylation annotation with nanopolish
      # If left blank, it will not be executed. And, with both parameters are set
      # it will automatically execute nanopolish to call methylation.
      # Only used with analysing a single genome
      # When analysing multiple genomes it must be set in the YAML file.
      # Check the example YAML at: https://github.com/fmalmeida/bacannot/blob/master/example_samplesheet.yaml

    --nanopolish_fast5 <string>                    Path to directory containing FAST5 files

    --nanopolish_fastq <string>                    Path to fastq files (file related to FAST5 files above)


""".stripIndent()
}
