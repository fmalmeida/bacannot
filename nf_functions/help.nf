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

                                  General Parameters

    --outdir <string>                              Output directory name

    --threads <int>                                Number of threads to use

    --bedtools_merge_distance                      By default, this process is not executed. For execution
                                                   one needs to provide a value.Minimum number of overlapping
                                                   bases for gene merge using bedtools merge. Negative values,
                                                   such as -20, means the number of required overlapping bases
                                                   for merging. Positive values, such as 5, means the maximum
                                                   distance accepted between features for merging.

                                  Configuring Input options
                (Users can give either a genome in FASTA file or raw reads in FASTQ)

    --genome <string>                              Set path to the genome in FASTA file. Users can annotate more than
                                                   one genome at once by using glob patters, such as "*.fasta"

                ( If used together at once, the different NGS reads, short and long reads,
                 must be from the same sample, one sample at a time. If you want to give
                 more than one sample at once you must use only one NGS read type as input
                 since we can't guaratee the order they are picked by nextflow )

    --sreads_paired                                Illumina paired end reads in fastq for assembly before annotation.

    --sreads_single                                Illumina unpaired reads in fastq for assembly before annotation.

    --lreads                                       Path to longreads in fastq assembly before annotation (ONT or Pacbio).

    --lreads_type                                  Tells the technology of the input longreads: [ nanopore or pacbio ].


                                  Prokka complementary parameters

    --prokka_kingdom <string>                      Prokka annotation mode. Possibilities (default 'Bacteria'):
                                                   Archaea|Bacteria|Mitochondria|Viruses

    --prokka_genetic_code <int>                    Genetic Translation code. Must be set if kingdom is not
                                                   default (in blank).

    --prokka_use_rnammer                           Tells prokka wheter to use rnammer instead of barrnap.


                                  Blast alignment parameters

    --blast_virulence_minid                        Min. identity % for virulence annotation. Default 90.

    --blast_virulence_mincov                       Min. gene coverage for virulence annotation. Default 80.

    --blast_resistance_minid                       Min. identity % for resistance annotation. Default 90.

    --blast_resistance_mincov                      Min. gene coverage for resistance annotation. Default 80.

    --blast_MGEs_minid                             Min. identity % for ICEs and prophage annotation. Default 65.

    --blast_MGEs_mincov                            Min. query coverage for ICEs and prophage annotation. Default 65.

    --plasmids_minid                               Min. identity % for plasmid detection. Default 90.

    --plasmids_mincov                              Min. query coverage for plasmid detection. Default 60.


                                  Configure resfinder optional parameter

    --resfinder_species                            It sets the species to be used for Resfinder annotation. If blank,
                                                   it will not be executed. Must be identical (without the *) as written
                                                   in their webservice https://cge.cbs.dtu.dk/services/ResFinder/.

                                  Configure Optional processes

    --skip_virulence_search                        Tells wheter you do not want to execute virulence annotation

    --skip_resistance_search                       Tells wheter you do not want to execute resistance annotation

    --skip_iceberg_search                          Tells wheter you do not want to execute ICE annotation

    --skip_prophage_search                         Tells wheter you do not want to execute prophage annotation

    --skip_plasmid_search                          Tells wheter you do not want to execute plasmid detection

    --skip_kofamscan                               Tells wheter you do not want to execute KO annotation with kofamscan


                            Configure optional Methylation annotation with nanopolish

                    ( If left blank, it will not be executed. And, with both parameters are set
                      it will automatically execute nanopolish to call methylation. For using
                      these parameters, the pipeline must be used with one sample at a time
                      since we can't guaratee the order the files are picked by nextflow )

    --nanopolish_fast5_dir <string>                Path to directory containing FAST5 files

    --nanopolish_fastq_reads <string>              Path to fastq files (file related to FAST5 files above)


""".stripIndent()
}
