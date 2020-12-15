def exampleMessage() {
   log.info """
   Example Usages:

      Simple Klebsiella genome annotation using all pipeline's optional annotation processes
\$ nextflow run fmalmeida/bacannot --threads 3 --outdir kp25X --genome kp_ont.contigs.fasta --bedtools_merge_distance -20 --blast_virulence_minid 90 \
--blast_virulence_mincov 80 --blast_MGEs_minid 70 --blast_MGEs_mincov 60 --nanopolish_fast5 ./fast5_pass --nanopolish_fastq ./kp_ont.fastq \
--resfinder_species "Klebsiella"

     Various genomes at once
\$ nextflow run fmalmeida/bacannot --threads 3 --outdir teste --genome "input/*.fasta" --resfinder_species "Escherichia coli"

     Various paired end shortreads samples at once
\$ nextflow run fmalmeida/bacannot --threads 3 --outdir teste --sreads_paired "*{1,2}.fastq" --resfinder_species "Klebsiella"

     Combining different NGS reads -- Must be used with only one sample at a time.
\$ nextflow run fmalmeida/bacannot --threads 3 --outdir teste_one_sample --sreads_paired "sample1_{1,2}.fastq" --lreads "sample1_nanopore.fastq" \
--lreads_type "nanopore" --sreads_single "sample1_sr_merged.fastq"


""".stripIndent()
}
