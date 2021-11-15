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

   Get template configuration file:

   nextflow run fmalmeida/bacannot [--get_config]

   Execution Reports:

   nextflow run fmalmeida/bacannot [OPTIONS] [-with-report] [-with-trace] [-with-timeline]

   OPTIONS:

      # Input configuration
      # Users can give either a genome in FASTA file or raw reads in FASTQ
      # The analysis is configured via a samplesheet
      # Check the example samplesheet at: https://github.com/fmalmeida/bacannot/blob/master/example_samplesheet.yaml
      #
      # Also documented at: https://bacannot.readthedocs.io/en/latest/samplesheet.html

    --input <string>                               Set path to the input samplesheet.

      # Annotation configuration
      # Read it and configure it properly

      # General Parameters

    --output <string>                              Output directory name

    --threads <int>                                Number of threads to use

    --parallel_jobs <int>                          Number of jobs to run in parallel. Each job can consume up
                                                   to N threads (--threads). If nothing is provided, it let's nextflow automatically handle it, which is the default.

    --bedtools_merge_distance <int>                By default, this process is not executed. For execution
                                                   one needs to provide a minimum number of overlapping
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

    --blast_virulence_minid <int>                  Min. identity % for virulence annotation. Default 90.

    --blast_virulence_mincov <int>                 Min. gene/subject coverage for virulence annotation. Default 80.

    --blast_resistance_minid <int>                 Min. identity % for resistance annotation. Default 90.

    --blast_resistance_mincov <int>                Min. gene/subject coverage for resistance annotation. Default 80.

    --blast_MGEs_minid <int>                       Min. identity % for ICEs and prophage annotation. Default 65.

    --blast_MGEs_mincov <int>                      Min. gene/subject coverage for ICEs and prophage annotation. Default 65.

    --plasmids_minid <int>                         Min. identity % for plasmid detection. Default 90.

    --plasmids_mincov <int>                        Min. coverage for plasmid detection. Default 60.

    --blast_custom_minid <int>                     Min. identity % for the annotation using user's custom database. Default 0.

    --blast_custom_mincov <int>                    Min. gene/subject coverage % for the annotation using user's custom database. Default 0.

      # User's custom database for annotation
      # Must be in gene nucleotide FASTA
      #
      # Well documented at: https://bacannot.readthedocs.io/en/latest/custom-db.html

    --custom_db <string>                           Path to the nucleotide FASTA file containing the user's custom database for annotation.
                                                   Multiple FASTAs can be provided separated by comma. E.g. db1.fasta,db2.fasta,...


      # Configure a default resfinder species panel for all samples
      # If a sample has another value inside the samplesheet, the pipeline will use 
      # the one found inside the samplesheet for that specific sample.
      #
      # Also documented at: https://bacannot.readthedocs.io/en/latest/samplesheet.html

    --resfinder_species <string>                   It sets the species to be used for Resfinder annotation. If blank,
                                                   it will not be executed. Must be identical (without the *) as written
                                                   in their webservice https://cge.cbs.dtu.dk/services/ResFinder/.
                                                   If your species is not available in Resfinder panels, you may use it 
                                                   with the "Other" panel (--resfinder_species "Other").
                                                   E.g. 'Escherichia coli'; 'Klebsiella' ...

      # Configure (on/off) optional processes

    --skip_virulence_search                        Tells whether you do not want to execute virulence annotation

    --skip_resistance_search                       Tells whether you do not want to execute resistance annotation

    --skip_iceberg_search                          Tells whether you do not want to execute ICE annotation

    --skip_prophage_search                         Tells whether you do not want to execute prophage annotation

    --skip_plasmid_search                          Tells whether you do not want to execute plasmid detection

    --skip_kofamscan                               Tells whether you do not want to execute KO annotation with kofamscan

""".stripIndent()
}
