def exampleMessage() {
   log.info """
   Example Usages:

      ## Launching interactive graphical interface
      nf-core launch fmalmeida/bacannot

      ## Single genome annotation, with customization
\$ nextflow run fmalmeida/bacannot --outdir TESTE --threads 3 --genome assembly.fasta --bedtools_merge_distance -20 --blast_virulence_minid 90 \
--blast_virulence_mincov 80 --blast_MGEs_minid 70 --blast_MGEs_mincov 60 --nanopolish_fast5 ./fast5_pass --nanopolish_fastq ./ont.fastq \
--resfinder_species "Klebsiella"

     ## Multiple genome annotation, with custom database annotation
     ## using either raw reads or assembled genomes
\$ nextflow run fmalmeida/bacannot --outdir TESTE --threads 3 --in_yaml samplesheet.yaml --custom_db db1.fasta

     ## Annotating from raw reads
\$ nextflow run fmalmeida/bacannot --sreads_paired "sample1_{1,2}.fastq" --lreads "sample1_lreads.fastq" --lreads_type nanopore \
--outdir TESTE --skip_kofamscan --threads 5 --nanopolish_fastq_reads "sample1_lreads.fastq" --nanopolish_fast5_dir "fast5_pass_dir"

    ## Running with a configuration file
\$ nextflow run fmalmeida/bacannot -c bacannot.config

""".stripIndent()
}
