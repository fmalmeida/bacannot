#!/usr/bin/env nextflow

/*
          Generic Pipeline for Prokariotic Genome Annotation

          This pipeline was built to make it easier for those
          with minimum bioinformatics expertise to perform a
          comprehensive genome annotation analysis.

          Feel free to add more processes to this code if need.

*/

def helpMessage() {
   log.info """
   Usage:
   nextflow run fmalmeida/bacannot [--help] [ -c nextflow.config ] [OPTIONS] [-with-report] [-with-trace] [-with-timeline]
   Comments:

   This pipeline contains a massive amount of configuration variables and its usage as CLI parameters would
   cause the command to be huge.
   Therefore, it is extremely recommended to use the nextflow.config configuration file in order to make
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

                            General Parameters - Mandatory

    --outDir <string>                              Output directory name
    --threads <int>                                Number of threads to use
    --genome <string>                              Query Genome file
    --bedtools_merge_distance                      Minimum number of overlapping bases for gene merge
                                                   using bedtools merge (default: 0)

                            Prokka complementary parameters

    --prokka_center <string>                       Your institude acronym to be used by prokka when
                                                   renaming contigs.
    --prokka_kingdom <string>                      Prokka annotation mode. Possibilities (default 'Bacteria'):
                                                   Archaea|Bacteria|Mitochondria|Viruses
    --prokka_genetic_code <int>                    Genetic Translation code. Must be set if kingdom is not
                                                   default (in blank).
    --prokka_use_rnammer                           Tells prokka wheter to use rnammer instead of barrnap.
    --prokka_genus <string>                        Set only if you want to search only a specific genus database

                            Diamond (blastx) search parameters

    --diamond_virulence_identity                   Min. identity % for virulence annotation
    --diamond_virulence_queryCoverage              Min. query coverage for virulence annotation
    --diamond_MGEs_identity                        Min. identity % for ICEs and prophage annotation
    --diamond_MGEs_queryCoverage                   Min. query coverage for ICEs and prophage annotation
    --diamond_minimum_alignment_length             Min. alignment length for diamond annotation

                            Configure Optional processes

    --not_run_virulence_search                     Tells wheter you want or not to execute virulence annotation
    --not_run_vfdb_search                          Tells wheter you want or not to used VFDB database for virulence
                                                   annotation. It is useless if virulence_search is not true
    --not_run_victors_search                       Tells wheter you want or not to used victors database for virulence
                                                   annotation. It is useless if virulence_search is not true
    --not_run_resistance_search                    Tells wheter you want or not to execute resistance annotation
    --not_run_iceberg_search                       Tells wheter you want or not to execute ICE annotation
    --not_run_prophage_search                      Tells wheter you want or not to execute prophage annotation
    --not_run_kofamscan                            Tells wheter you want or not to execute KO annotation with kofamscan

                            Configure Optional Pangenome analysis with Roary

    --roary_reference_genomes <string>             Used to set path to reference genomes to be used in the pangenome
                                                   analysis with Roary. Whenever set, the pipeline will automatically
                                                   execute Roary pangenome analysis. Example: "path/reference/*.fasta"
                                                   They must be all in one directory and they must no be links. They
                                                   must be the hard file.

                      Configure optional Methylation annotation with nanopolish
                      If left blank, it will not be executed. And, with both parameters are set
                      it will automatically execute nanopolish to call methylation

    --nanopolish_fast5_dir <string>                Path to directory containing FAST5 files
    --nanopolish_fastq_reads <string>              Path to fastq files (file related to FAST5 files above)


""".stripIndent()
}

def exampleMessage() {
   log.info """
   Example Usages:
      Simple Klebsiella genome annotation using all pipeline's optional annotation processes
nextflow run fmalmeida/bacannot --threads 3 --outDir kp25X --genome Kp31_BC08.contigs.fasta --bedtools_merge_distance -20 --prokka_center UNB --diamond_virulence_identity 90 --diamond_virulence_queryCoverage 80 --diamond_MGEs_identity 70 --diamond_MGEs_queryCoverage 60 --diamond_minimum_alignment_length 200 --virulence_search --vfdb_search --victors_search --resistance_search --ice_search --prophage_search --execute_kofamscan --nanopolish_fast5_dir fast5_pass --nanopolish_fastq_reads Kp31_BC08.fastq


""".stripIndent()
}

/*
                                                  Display Help Message
*/

params.help = false
 // Show help emssage
 if (params.help){
   helpMessage()
   //file('work').deleteDir()
   //file('.nextflow').deleteDir()
   exit 0
}

// CLI examples
params.examples = false
 // Show help emssage
 if (params.examples){
   exampleMessage()
   exit 0
}

// Get configuration file
params.get_config = false
if (params.get_config) {
  new File("bacannot.config") << new URL ("https://github.com/fmalmeida/bacannot/raw/master/nextflow.config").getText()
  println ""
  println "bacannot.config file saved in working directory"
  println "After configuration, run:"
  println "nextflow run fmalmeida/bacannot -c ./bacannot.config"
  println "Nice code!\n"

  //file('work').deleteDir()
  //file('.nextflow').deleteDir()
  exit 0
}

/*

                                  Setting default parameters

*/

params.prefix = 'out'
params.outDir = 'outdir'
params.threads = 2
params.bedtools_merge_distance = 0
params.prokka_center = 'Centre'
params.prokka_kingdom = ''
params.prokka_genetic_code = false
params.prokka_use_rnammer = false
params.prokka_genus = ''
params.diamond_virulence_identity = 90
params.diamond_virulence_queryCoverage = 90
params.diamond_MGEs_identity = 65
params.diamond_MGEs_queryCoverage = 65
params.diamond_minimum_alignment_length = 200
params.not_run_virulence_search = false
params.not_run_vfdb_search = false
params.not_run_victors_search = false
params.not_run_resistance_search = false
params.not_run_iceberg_search = false
params.not_run_prophage_search = false
params.not_run_kofamscan = false
params.roary_reference_genomes = false

/*

                                  Loading general parameters

*/

genome = file(params.genome)
reference_genomes = (params.roary_reference_genomes) ? Channel.fromPath( params.roary_reference_genomes ) : Channel.empty()
prefix = params.prefix
outDir = params.outDir
threads = params.threads
// This parameter sets the minimum number of overlapping bases for gene merge.
bedDistance = params.bedtools_merge_distance
// Diamond (blastx) parameters
diamond_virulence_identity = params.diamond_virulence_identity
diamond_virulence_queryCoverage = params.diamond_virulence_queryCoverage
diamond_MGEs_identity = params.diamond_MGEs_identity
diamond_MGEs_queryCoverage = params.diamond_MGEs_queryCoverage
diamond_minimum_alignment_length = params.diamond_minimum_alignment_length
params.virulence_search = (params.not_run_virulence_search) ? false : true
params.vfdb_search = (params.not_run_vfdb_search) ? false : true
params.victors_search = (params.not_run_victors_search) ? false : true
params.resistance_search = (params.not_run_resistance_search) ? false : true
params.iceberg_search = (params.not_run_iceberg_search) ? false : true
params.prophage_search = (params.not_run_prophage_search) ? false : true
params.execute_kofamscan = (params.not_run_kofamscan) ? false : true

/*

                              Pipeline execution begins

*/

process MLST {
   publishDir "${outDir}/MLST", mode: 'copy'
   container = 'fmalmeida/bacannot:latest'

   input:
   file input from genome

   output:
   file "${prefix}_mlst_analysis.txt" optional true
   file "${prefix}_novel_alleles.fasta" optional true

   script:
   """
   source activate MLST ;
   mlst --quiet --novel ${prefix}_novel_alleles.fasta $input > ${prefix}_mlst_analysis.txt
   """
}

process prokka {
    publishDir outDir, mode: 'copy'
    container = 'fmalmeida/bacannot:latest'

    input:
    file input from genome

    output:
    file "prokka/${prefix}.gff" into annotation_gff_prokka, annotation_gff_prokka_roary
    file "prokka/${prefix}.gbk" into annotation_gbk_prokka
    file "prokka/"
    file "prokka/${prefix}*.fna" into renamed_genome
    file "prokka/${prefix}*.faa" into genes_aa_sql, genes_aa_kofamscan, amrfinder_input
    file "prokka/${prefix}*.ffn" into genes_sequences_vfdb, genes_sequences_iceberg, genes_nt_sql, genes_sequences_phast, genes_sequences_victors

    script:
    kingdom = (params.prokka_kingdom) ? "--kingdom ${params.prokka_kingdom}" : ''
    gcode = (params.prokka_genetic_code) ? "--gcode ${params.prokka_genetic_code}" : ''
    rnammer = (params.prokka_use_rnammer) ? "--rnammer" : ''
    genus = (params.prokka_genus) ? "--genus ${params.prokka_genus} --usegenus" : ''
    """
    source activate PROKKA ;
    prokka $kingdom $gcode $rnammer --outdir prokka --cpus $threads --centre ${params.prokka_center} \
    --mincontiglen 200 $genus --prefix $prefix $input
    """
}

process create_roary_input {
  container = 'fmalmeida/bacannot:latest'

  input:
  file references from reference_genomes

  output:
  file "${references.baseName}/${references.baseName}.gff" into roary_inputs

  when:
  (params.reference_genomes)

  script:
  """
  ## First, run prokka on all genomes
  source activate PROKKA ;
  prokka --cpus $threads --centre ${params.prokka_center} \
  --mincontiglen 200 --prefix ${references.baseName} $references ;
  """
}

process roary_pangenome {
  publishDir "${outDir}/Roary_pangenome", mode: 'copy'
  container = 'fmalmeida/bacannot:latest'

  input:
  file gffs from roary_inputs.mix(annotation_gff_prokka_roary).collect()

  output:
  file "roary_output"
  file "roary_inputs"

  when:
  (params.reference_genomes)

  script:
  """
  ## Activate environment
  source activate ROARY ;

  ## Inputs
  mkdir roary_inputs ;
  cp ${gffs} roary_inputs ;

  ## Run roary pipeline
  roary -e -f roary_output --mafft -p $threads $gffs ;
  FastTree -nt -gtr roary_output/core_gene_alignment.aln > roary_output/output.newick

  ## Entering Dir
  cd roary_output ;

  ## Plotting
  roary2svg.pl gene_presence_absence.csv > pan_genome.svg ;
  conda deactivate ;
  source activate ROARY_PLOTS ;
  python /usr/local/bin/roary_plots.py output.newick gene_presence_absence.csv ;
  """
}

process rRNA {
   publishDir "${outDir}/rRNA", mode: 'copy'
   container = 'fmalmeida/bacannot:latest'

   input:
   file input from renamed_genome

   output:
   file "${prefix}_rRNA.gff" into rrna_gff
   file "${prefix}_rRNA.fa" optional true

   script:
   """
   barrnap -o ${prefix}_rRNA.fa < $input > ${prefix}_rRNA.gff
   """
}

/*

                        Remove genome sequence and comments from GFF.
                        Also, it masks the genome file.

*/

process masking_genome {
  container 'fmalmeida/bacannot:latest'

  input:
  file input from renamed_genome
  file 'gff' from annotation_gff_prokka

  output:
  file "${prefix}_clear.gff" into clear_gff
  file "${prefix}_masked_genome.fasta" into masked_genome_vfdb, masked_genome_phast

  """
  grep "ID=" gff | awk '{ print \$1 "\t" \$4 "\t" \$5 }' > cds_prokka.bed ;
  grep "ID=" gff > ${prefix}_clear.gff ;
  maskFastaFromBed -fi $input -fo ${prefix}_masked_genome.fasta -bed cds_prokka.bed
  """
}

/*

                                                  Compute GC content to plot in JBrowse

*/

process compute_GC {
  container 'fmalmeida/bacannot:latest'

  input:
  file 'input.fasta' from renamed_genome

  output:
  file "input_GC_500_bps.sorted.bedGraph" into gc_content_jbrowse
  file "input.sizes" into gc_sizes_jbrowse

  """
  # Index
  samtools faidx input.fasta ;
  # Take Sizes
  cut -f 1,2 input.fasta.fai > input.sizes ;
  # Create sliding window
  bedtools makewindows -g input.sizes -w 500 > input_500_bps.bed ;
  # Compute GC
  bedtools nuc -fi input.fasta -bed input_500_bps.bed > input_500_bps_nuc.txt ;
  # Create bedGraph
  awk 'BEGIN{FS="\\t"; OFS="\\t"} FNR > 1 { print \$1,\$2,\$3,\$5 }' input_500_bps_nuc.txt > input_GC_500_bps.bedGraph
  # Sort
  bedtools sort -i input_GC_500_bps.bedGraph > input_GC_500_bps.sorted.bedGraph
  """
}

/*

                KOfamscan process. Here we search our predicted proteins in order to retrieve KO (KEGG)
                Information about then. This makes easier to plot Methabolic Fluxes in KEGG.

*/

process kofamscan {
  publishDir outDir, mode: 'copy'
  container = 'fmalmeida/bacannot:kofamscan'
  x = "Executing KOfamscan - Its output can be directly put in KEGG Mapper for visualization"
  tag { x }

  input:
  file 'proteins.faa' from genes_aa_kofamscan

  output:
  file "KOfamscan/${prefix}_ko_*"
  file "KOfamscan/${prefix}_ko_forKEGGMapper.txt" into kofamscan_hits

  when:
  ( params.execute_kofamscan )

  script:
  """
  mkdir KOfamscan ;
  kofamscan -o KOfamscan/${prefix}_ko_detailed.txt --cpu=${threads} proteins.faa ;
  kofamscan -o KOfamscan/${prefix}_ko_forKEGGMapper.txt --cpu=${threads} -f mapper-one-line proteins.faa ;
  """
}

/*

                      Initiaing DIAMOND searches against specific databases:
                                VFDB, Victors,ICEberg and Phast.

*/

process vfdb {
  if ( params.virulence_search && params.vfdb_search ) {
  publishDir outDir, mode: 'copy',
  saveAs: {filename ->
  //This line saves the files with specific sufixes in specific folders
  if (filename.indexOf(".tsv") > 0 ) "blasts/$filename"
  else if (filename.indexOf(".gff") > 0 ) "gffs/only_against_maskedGenome/$filename"
  }}
  container 'fmalmeida/bacannot:latest'
  x = ( params.virulence_search && params.vfdb_search
        ? "Process is being executed"
        : "Process was skipped by the user")
  tag { x }

  input:
  file genome from masked_genome_vfdb
  file genes from genes_sequences_vfdb

  output:
  file "${prefix}_vfdb.gff" into annotation_gff_vfdb
  file "*.tsv"
  file "virulence_vfdb_predictedGenes.tsv" into vfdb_blast_genes, vfdb_blast_genes2

  script:
  if ( params.virulence_search && params.vfdb_search )
  """
  # First step, with masked genome
  diamond blastx --query-gencode 11 --db /work/vfdb/vfdb_prot -o blast_result.tmp \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore stitle\
  --query $genome --query-cover $diamond_virulence_queryCoverage ;

  ## Convert it to gff
  awk -v len=$diamond_minimum_alignment_length '{if (\$4 >= len) print }' blast_result.tmp > virulence_VFDB_maskedGenome.tsv ;
  python2 /usr/local/bin/blast2gff.py -b virulence_VFDB_maskedGenome.tsv -p vfdb -t virulence -F > ${prefix}_vfdb.gff ;

  # Second step, with predicted genes
  diamond blastx --query-gencode 11 --db /work/vfdb/vfdb_prot -o blast_result_genes.tmp \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore stitle \
  --query $genes --query-cover $diamond_virulence_queryCoverage ;

  ## Convert it to gff
  awk -v id=$diamond_virulence_queryCoverage '{if (\$3 >= id) print }' blast_result_genes.tmp > virulence_vfdb_predictedGenes.tsv
  """
  else
  """
  touch virulence_vfdb_predictedGenes.tsv ;
  touch ${prefix}_vfdb.gff
  """

}

process victors {
  if ( params.virulence_search && params.victors_search ) {
  publishDir outDir, mode: 'copy',
  saveAs: {filename ->
  //This line saves the files with specific sufixes in specific folders
  if (filename.indexOf(".tsv") > 0 ) "blasts/$filename"
  else if (filename.indexOf(".gff") > 0 ) "gffs/only_against_maskedGenome/$filename"
}}
  container 'fmalmeida/bacannot:latest'
  x = ( params.virulence_search && params.victors_search
        ? "Process is being executed"
        : "Process was skipped by the user")
  tag { x }

  input:
  file genes from genes_sequences_victors

  output:
  file "virulence_victors_predictedGenes.tsv" into victors_blast_genes, victors_blast_genes2

  script:
  if ( params.virulence_search && params.victors_search )
  """
  # Blast predicted genes
  diamond blastx --query-gencode 11 --db /work/victors/victors_prot -o blast_result_genes.tmp \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore stitle \
  --query $genes --query-cover $diamond_virulence_queryCoverage ;

  ## Convert it to gff
  awk -v id=$diamond_virulence_queryCoverage '{if (\$3 >= id) print }' blast_result_genes.tmp > virulence_victors_predictedGenes.tsv
  """
  else
  """
  touch virulence_victors_predictedGenes.tsv
  """
}

process amrfinder {
  if ( params.resistance_search ) {
  publishDir "${outDir}/resistance/AMRFinderPlus_table", mode: 'copy' }
  container 'fmalmeida/bacannot:latest'
  x = ( params.resistance_search
        ? "Process is being executed"
        : "Process was skipped by the user")
  tag { x }

  input:
  file proteins from amrfinder_input

  output:
  file "AMRFinder_output.tsv" into amrfinder_output optional true

  when:
  ( params.resistance_search )

  script:
  """
  source activate AMRFINDERPLUS ;
  amrfinder -p $proteins --plus -o AMRFinder_output.tsv
  """
}

process rgi_annotation {
  if ( params.resistance_search ) {
  publishDir "${outDir}/resistance/RGI_annotation", mode: 'copy' }
  container 'fmalmeida/bacannot:latest'
  x = ( params.resistance_search
        ? "Process is being executed"
        : "Process was skipped by the user")
  tag { x }

  input:
  file input from renamed_genome

  output:
  file "Perfect_RGI_${params.prefix}_hits.txt" into rgi_perfect optional true
  file "Strict_RGI_${params.prefix}_hits.txt" into rgi_strict optional true
  file "RGI_${params.prefix}.txt" into rgi_output optional true
  file "*RGI_${params.prefix}*" optional true

  when:
  ( params.resistance_search )

  script:
  """
  source activate RGI ;
  rgi main -i $input -o RGI_${params.prefix} -t contig -a DIAMOND -n ${params.threads} --exclude_nudge --clean ;

  ## Parse perfect hits
  sed 's/ # /#/g' RGI_${params.prefix}.txt \
  | awk 'BEGIN { FS = "\t"; OFS="\\t" } ; { print \$2,\$3,\$4,\$5,\$6,\$9,\$11,\$15,\$16,\$17 }' \
  | awk '{ if (\$5 == "Perfect") print }' > Perfect_RGI_${params.prefix}_hits.txt

  ## Parse strict hits
  sed 's/ # /#/g' RGI_${params.prefix}.txt \
  | awk 'BEGIN { FS = "\t"; OFS="\\t" } ; { print \$2,\$3,\$4,\$5,\$6,\$9,\$11,\$15,\$16,\$17 }' \
  | awk '{ if (\$5 == "Strict") print }' > Strict_RGI_${params.prefix}_hits.txt
  """
}

process phast {
  if ( params.prophage_search ) {
  publishDir outDir, mode: 'copy',
  saveAs: {filename ->
  //This line saves the files with specific sufixes in specific folders
  if (filename.indexOf(".tsv") > 0 ) "blasts/$filename"
  else if (filename.indexOf(".gff") > 0 ) "gffs/only_against_maskedGenome/$filename"
}}
  container 'fmalmeida/bacannot:latest'
  x = ( params.prophage_search
        ? "Process is being executed"
        : "Process was skipped by the user")
  tag { x }

  input:
  file genome from masked_genome_phast
  file genes from genes_sequences_phast

  output:
  file "${prefix}_phast.gff" into annotation_gff_phast
  file "*.tsv"
  file "prophage_phast_predictedGenes.tsv" into phast_blast_genes, phast_blast_genes2

  script:
  if ( params.prophage_search )
  """
  # First step, with masked genome
  diamond blastx --query-gencode 11 --db /work/phast/phast_prot -o blast_result.tmp \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore stitle\
  --query $genome --query-cover $diamond_MGEs_queryCoverage ;

  ## Convert it to gff
  awk -v len=$diamond_minimum_alignment_length '{ if (\$4 >= len) print }' blast_result.tmp > prophage_phast_maskedGenome.tsv ;
  python2 /usr/local/bin/blast2gff.py -b prophage_phast_maskedGenome.tsv -p phast -t prophage -F > ${prefix}_phast.gff ;

  # Second step, with predicted genes
  diamond blastx --query-gencode 11 --db /work/phast/phast_prot -o blast_result_genes.tmp \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore stitle \
  --query $genes --query-cover $diamond_MGEs_queryCoverage ;

  ## Convert it to gff
  awk -v id=$diamond_MGEs_identity '{if (\$3 >= id) print }' blast_result_genes.tmp > prophage_phast_predictedGenes.tsv
  """
  else
  """
  touch prophage_phast_predictedGenes.tsv ;
  touch ${prefix}_phast.gff
  """
}

process phigaro {
  if ( params.prophage_search ) {
  publishDir "${outDir}/prophages/phigaro", mode: 'copy'}
  container 'fmalmeida/bacannot:latest'
  x = ( params.prophage_search
        ? "Process is being executed"
        : "Process was skipped by the user")
  tag { x }

  input:
  file "assembly.fasta" from renamed_genome

  output:
  file "assembly.phg" into phigaro_txt
  file "prophages.bed" into phigaro_bed
  file "assembly.phg.html"

  when:
  ( params.prophage_search )

  script:
  """
  touch assembly.phg assembly.phg.html ;
  seqtk seq -L 20000 assembly.fasta > assembly-L20000.fasta ;
  phigaro -f assembly-L20000.fasta -c /work/phigaro/config.yml -e html txt -o assembly.phg -p --not-open ;

  # Create BED
  grep -v "taxonomy" assembly.phg | awk 'BEGIN { FS = "\t"; OFS="\\t" } { print \$1,\$2,\$3 }' > prophages.bed
  """

}

process find_GIs {
  publishDir "${outDir}/predicted_GIs", mode: 'copy'
  container 'fmalmeida/bacannot:latest'

  input:
  file "annotation.gbk" from annotation_gbk_prokka

  output:
  file "predicted_GIs.bed" into predicted_GIs optional true

  script:
  """
  source activate find_GIs ;
  python /work/pythonScripts/splitgenbank.py annotation.gbk && rm annotation.gbk ;
  for file in \$(ls *.gbk); do grep -q "CDS" \$file && Dimob.pl \$file \${file%%.gbk}_GIs.txt 2> dimob.err ; done
  for GI in \$(ls *.txt); do \
    awk -v contig="\$( echo \"gnl|${params.prokka_center}|\${GI%%_GIs.txt}\" )" \
    'BEGIN { FS = "\t"; OFS="\\t" } { print contig,\$2,\$3 }' \$GI >> predicted_GIs.bed ; \
  done
  """
}

process iceberg {
  if ( params.iceberg_search ) {
  publishDir outDir, mode: 'copy',
  saveAs: {filename ->
  //This line saves the files with specific sufixes in specific folders
  if (filename.indexOf(".tsv") > 0 ) "blasts/$filename"
  else if (filename.indexOf(".gff") > 0 ) "gffs/only_against_maskedGenome/$filename" }}
  container 'fmalmeida/bacannot:latest'
  x = ( params.iceberg_search
        ? "Process is being executed"
        : "Process was skipped by the user")
  tag { x }

  input:
  file genes from genes_sequences_iceberg

  output:
  file "ice_iceberg_predictedGenes.tsv" into iceberg_blast_genes, ice_blast_genes2

  script:
  if ( params.iceberg_search )
  """
  # Blast search
  diamond blastx --query-gencode 11 --db /work/iceberg/iceberg_prot -o blastprot2_result.tmp \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore stitle \
  --query $genes --query-cover $diamond_MGEs_queryCoverage ;

  ## Group it all
  awk -v id=$diamond_MGEs_identity '{if (\$3 >= id) print }' blastprot2_result.tmp > ice_iceberg_predictedGenes.tsv
  """
  else
  """
  touch ice_iceberg_predictedGenes.tsv ;
  """
}

/*
 * This process adds to the GFF the information about the additional
 * databases that were BLAST searched (with predicted genes). Note,
 * if the user has set its BLAST to FALSE, the blank file will be ignored
 * by the following process and will only perform tasks with those that
 * some content in it.
 */

process genes_blasted_to_gff {
  publishDir "${outDir}/gffs/only_against_predicted_genes", mode: 'copy'
  container 'fmalmeida/bacannot:renv'

  input:
  file 'gff' from clear_gff
  file 'RGI_output.txt' from rgi_output.ifEmpty('RGI_empty')
  file blastVFDB from vfdb_blast_genes
  file blastVictors from victors_blast_genes
  file blastIce from iceberg_blast_genes
  file blastPhast from phast_blast_genes
  file "AMRFinder_output.tsv" from amrfinder_output.ifEmpty('AMRFinder_empty')
  file 'kofamscan.txt' from kofamscan_hits.ifEmpty('kofamscan_empty')

  output:
  file "${prefix}_blast_genes.gff" into blast_genes_gff

  script:
  """
  ## Reformat KOfamscan Output
  while read line ; do id=\$(echo \$line | awk '{print \$1}') ; ko=\$(echo \$line | awk '{\$1=""; print \$0}' | \
  sed 's/\\s//' | sed 's/\\s/,/g') ; echo -e "\$id\\t\$ko" ; done < kofamscan.txt > formated.txt ;

  ## Add features
  addBlast2Gff.R -i $blastVFDB -g gff -o gff -d vfdb -t virulence -c ${params.diamond_virulence_queryCoverage};
  [ ! -s kofamscan.txt ] || addKO2Gff.R -i formated.txt -g gff -o gff -d KEGG ;
  addBlast2Gff.R -i $blastVictors -g gff -o gff -d victors -t virulence -c ${params.diamond_virulence_queryCoverage};
  addBlast2Gff.R -i $blastIce -g gff -o gff -d iceberg -t ice -c ${params.diamond_MGEs_queryCoverage};
  addBlast2Gff.R -i $blastPhast -g gff -o gff -d phast -t prophage -c ${params.diamond_MGEs_queryCoverage};
  [ ! -s RGI_output.txt ] || addRGI2gff.R -g gff -i RGI_output.txt -o gff ;
  [ ! -s AMRFinder_output.tsv ] || addNCBIamr2Gff.R -g gff -i AMRFinder_output.tsv -o ${prefix}_blast_genes.gff -t resistance -d AMRFinderPlus ;
  [ -s AMRFinder_output.tsv ] || mv gff ${prefix}_blast_genes.gff ;
  """
}

/*
 * This process will only merge gff files that have content in it. So,
 * gff files in blank from addional databases that were not searched
 * (with masked genome) will be ignored.
 */

process merge_gffs {
  publishDir "${outDir}/gffs/merged", mode: 'copy'
  container 'fmalmeida/bacannot:latest'

  input:
  file 'prokka_gff' from blast_genes_gff
  file 'vfdb_gff' from annotation_gff_vfdb
  file 'phast_gff' from annotation_gff_phast

  output:
  file "${prefix}_merged.gff" into merged_gff

  """
  # GFF from prokka
  [ ! -s prokka_gff ] || sed 's/ /_/g' prokka_gff | bedtools sort \
  | bedtools merge -d $bedDistance -s -c 2,3,6,7,8,9 -o distinct,distinct,max,distinct,distinct,distinct \
  | awk 'BEGIN { FS = "\t"; OFS="\\t" } { print \$1,\$4,\$5,\$2+1,\$3,\$6,\$7,\$8,\$9}' >> total_gff ;

  # GFF from VFDB
  [ ! -s vfdb_gff ] || sed 's/ /_/g' vfdb_gff | bedtools sort \
  | bedtools merge -d $bedDistance -s -c 2,3,6,7,8,9 -o distinct,distinct,max,distinct,distinct,distinct \
  | awk 'BEGIN { FS = "\t"; OFS="\\t" } { print \$1,\$4,\$5,\$2+1,\$3,\$6,\$7,\$8,\$9}' >> total_gff ;

  # GFF from PHAST
  [ ! -s phast_gff ] || sed 's/ /_/g' phast_gff | bedtools sort \
  | bedtools merge -d $bedDistance -s -c 2,3,6,7,8,9 -o distinct,distinct,max,distinct,distinct,distinct \
  | awk 'BEGIN { FS = "\t"; OFS="\\t" } { print \$1,\$4,\$5,\$2+1,\$3,\$6,\$7,\$8,\$9}' >> total_gff ;

  # Add header to it
  echo \"##gff-version 3\" > ${prefix}_merged.gff ;
  bedtools sort -i total_gff | bedtools merge -d $bedDistance -s -c 2,3,6,7,8,9 -o distinct,distinct,max,distinct,distinct,distinct \
  | awk 'BEGIN { FS = "\t"; OFS="\\t" } sub(",",";",\$9) { print \$1,\$4,\$5,\$2+1,\$3,\$6,\$7,\$8,\$9}' >> ${prefix}_merged.gff
  """
}

/*
 * Here, it starts the execution of summarization processes.
 */

process write_summary_tables {
  publishDir outDir, mode: 'copy',
  saveAs: {filename ->
  //This line saves the files with specific sufixes in specific folders
  if (filename.indexOf(".tsv") > 0 ) "/tmp/$filename"
  else if (filename.indexOf(".gff") > 0 ) "gffs/final/$filename" }
  container 'fmalmeida/bacannot:renv'

  input:
  file gff from merged_gff

  output:
  file "${prefix}_final.gff" into final_gff
  file "${prefix}_virulence.gff" into virulence_gff optional true
  file "${prefix}_resistance.gff" into resistance_gff optional true
  file "${prefix}_ices.gff" into ices_gff optional true
  file "${prefix}_prophage.gff" into prophage_gff optional true
  file "${prefix}_conjugation.gff" into conjugation_gff optional true
  file "${prefix}_efflux.gff" into efflux_gff optional true
  file "${prefix}_noHypothetical.gff" into noHypothetical_gff optional true
  file "${prefix}_hypothetical.gff" into hypothetical_gff optional true
  file "${prefix}_transposase.gff" into transposase_gff optional true
  file "${prefix}_amrfinderplus.gff" into amrfinderplus_gff optional true
  file "${prefix}_rgi.gff" into rgi_gff optional true
  file "${prefix}_amrfinderplus.tsv" into amrfinderplus_table optional true

  """
  # Reduce repeated values
  reduceRepeatedValues.R -i $gff -o tmp.gff ;
  sed -e 's/\\.,0/0/g' -e 's/0,\\./0/g' tmp.gff > gff ;
  tolower.R -i gff -o ${prefix}_final.gff;

  # Create gffs based on pattern
  touch ${prefix}_virulence.gff ${prefix}_resistance.gff ${prefix}_prophage.gff ${prefix}_ices.gff ;
  [[ \$(grep -c "virulence" ${prefix}_final.gff) -eq 0 ]] || awk '/virulence/ || /victors/ || vfdb' ${prefix}_final.gff > ${prefix}_virulence.gff ;
  [[ \$(grep -c "resistance" ${prefix}_final.gff) -eq 0 ]] || awk '/resistance/' ${prefix}_final.gff > ${prefix}_resistance.gff ;
  [[ \$(grep -c "AMRFinderPlus" ${prefix}_final.gff) -eq 0 ]] || awk '/AMRFinderPlus/' ${prefix}_final.gff > ${prefix}_amrfinderplus.gff ;
  [[ \$(grep -c "RGI" ${prefix}_final.gff) -eq 0 ]] || awk '/CARD-RGI/' ${prefix}_final.gff > ${prefix}_rgi.gff ;
  [[ \$(grep -c "prophage" ${prefix}_final.gff) -eq 0 ]] || awk '/prophage/ || \$2 ~ /phast/' ${prefix}_final.gff > ${prefix}_prophage.gff ;
  [[ \$(grep -c "ice" ${prefix}_final.gff) -eq 0 ]] || grep "iceberg" ${prefix}_final.gff > ${prefix}_ices.gff ;
  [[ \$(grep -c "efflux" ${prefix}_final.gff) -eq 0 ]] || grep "efflux" ${prefix}_final.gff > ${prefix}_efflux.gff ;
  [[ \$(grep -c -E "conjugative|conjugation|pilus|relaxase" ${prefix}_final.gff) -eq 0 ]] || grep -E "conjugative|conjugation|pilus|relaxase" ${prefix}_final.gff > ${prefix}_conjugation.gff ;
  [[ \$(grep -c "hypothetical" ${prefix}_final.gff) -eq 0 ]] || grep -v "hypothetical" ${prefix}_final.gff > ${prefix}_noHypothetical.gff ;
  [[ \$(grep -c "hypothetical" ${prefix}_final.gff) -eq 0 ]] || grep "hypothetical" ${prefix}_final.gff > ${prefix}_hypothetical.gff ;
  [[ \$(grep -c "transposase" ${prefix}_final.gff) -eq 0 ]] || grep "transposase" ${prefix}_final.gff > ${prefix}_transposase.gff ;

  # Write summary tables
  [ ! -s ${prefix}_amrfinderplus.gff ] || write_table_from_gff.R -i ${prefix}_amrfinderplus.gff -o ${prefix} -t amrfinderplus &> /tmp/error.txt ;
  """
}

process gff_to_gbk {
  publishDir "${outDir}/genbankFile", mode: 'copy'
  container 'fmalmeida/bacannot:latest'

  input:
  file gff from final_gff
  file input from genome

  output:
  file "*.genbank"

  """
  seqret -sequence $input -feature -fformat gff -fopenfile $gff -osformat genbank \
  -osname_outseq ${prefix} -ofdirectory_outseq gbk_file -auto
  """
}

/*
 * Calling Methylation with Nanopolish Call-Methylation
 */

//Get Files Necessary For Methylation Calling
if (params.fast5_dir && params.fastq_reads) {
  fast5 = Channel.fromPath( params.fast5_dir )
  nanopolish_lreads = file(params.fastq_reads)
  fast5_dir = Channel.fromPath( params.fast5_dir, type: 'dir' )
} else {
  fast5 = ''
  nanopolish_lreads = ''
  fast5_dir = ''
  }

process call_methylation {
  if (params.fast5_dir && params.fastq_reads) {
  publishDir "${outDir}/methylation", mode: 'copy' }
  container 'fmalmeida/bacannot:latest'
  x = ( params.fast5_dir && params.fastq_reads
        ? "Methylated sites are being calculated"
        : "Process was skipped by the user")
  tag { x }

  input:
  file 'input.fa' from renamed_genome
  file 'reads.fq' from nanopolish_lreads
  file fast5
  val fast5_dir from fast5_dir

  output:
  file "*_calls.tsv" optional true
  file "*_frequency.tsv" optional true
  file "cpg_frequency.bedGraph" into cpg_bedGraph
  file "gpc_frequency.bedGraph" into gpc_bedGraph
  file "dam_frequency.bedGraph" into dam_bedGraph
  file "dcm_frequency.bedGraph" into dcm_bedGraph
  file "chr.sizes"  into chr_sizes

  script:
  if (params.fast5_dir && params.fastq_reads)
  """
  # Index Our Fast5 Data
  nanopolish index -d "${fast5_dir}" reads.fq ;

  # Map Our Indexed Reads to Our Genome
  minimap2 -a -x map-ont input.fa reads.fq | samtools sort -T tmp -o reads_output.sorted.bam ;
  samtools index reads_output.sorted.bam ;

  # Call Methylation
  ## cpg , gpc , dam and dcm
  nanopolish call-methylation -r reads.fq -b reads_output.sorted.bam -g input.fa --methylation cpg > cpg_calls.tsv ;
	nanopolish call-methylation -r reads.fq -b reads_output.sorted.bam -g input.fa --methylation gpc > gpc_calls.tsv ;
	nanopolish call-methylation -r reads.fq -b reads_output.sorted.bam -g input.fa --methylation dam > dam_calls.tsv ;
	nanopolish call-methylation -r reads.fq -b reads_output.sorted.bam -g input.fa --methylation dcm > dcm_calls.tsv ;

  # Calculate Methylation Frequencies
  ## cpg , gpc , dam and dcm
  /work/nanopolish_scripts/calculate_methylation_frequency.py cpg_calls.tsv > cpg_frequency.tsv ;
  /work/nanopolish_scripts/calculate_methylation_frequency.py gpc_calls.tsv > gpc_frequency.tsv ;
  /work/nanopolish_scripts/calculate_methylation_frequency.py dam_calls.tsv > dam_frequency.tsv ;
  /work/nanopolish_scripts/calculate_methylation_frequency.py dcm_calls.tsv > dcm_frequency.tsv ;

  # Transform These TSV files into bedGraph
  ## cpg
  [ ! -s cpg_frequency.tsv ] || grep -v "start" cpg_frequency.tsv | \
  awk '{ print \$1 "\t" \$2 "\t" \$3 "\t" \$7 }' > cpg_frequency.bedGraph ;
  ## gpc
  [ ! -s gpc_frequency.tsv ] || grep -v "start" gpc_frequency.tsv | \
  awk '{ print \$1 "\t" \$2 "\t" \$3 "\t" \$7 }' > gpc_frequency.bedGraph ;
  ## dam
  [ ! -s dam_frequency.tsv ] || grep -v "start" dam_frequency.tsv |
  awk '{ print \$1 "\t" \$2 "\t" \$3 "\t" \$7 }' > dam_frequency.bedGraph ;
  ## dcm
  [ ! -s dcm_frequency.tsv ] || grep -v "start" dcm_frequency.tsv |
  awk '{ print \$1 "\t" \$2 "\t" \$3 "\t" \$7 }' > dcm_frequency.bedGraph ;

  # Create Contig Sizes File
  seqtk comp input.fa | awk '{ print \$1 "\t" \$2 }' > chr.sizes
  """
  else
  """
  touch cpg_frequency.bedGraph gpc_frequency.bedGraph dam_frequency.bedGraph dcm_frequency.bedGraph chr.sizes
  """
}

process jbrowse {
  publishDir "${outDir}/jbrowse/", mode: 'copy'
  container 'fmalmeida/bacannot:jbrowse'

  input:
  file input from renamed_genome
  file gff from final_gff
  file 'GC_content.bedGraph' from gc_content_jbrowse
  file 'GC_content.sizes' from gc_sizes_jbrowse
  file 'rrna.gff' from rrna_gff
  file 'resistance' from resistance_gff
  file 'amrfinder' from amrfinderplus_gff
  file 'rgi' from rgi_gff
  file 'virulence' from virulence_gff
  file 'prophage' from prophage_gff
  file 'prophages.bed' from phigaro_bed
  file 'all_GIs.bed' from predicted_GIs
  file 'ices' from ices_gff
  file 'conjugation' from conjugation_gff
  file 'efflux' from efflux_gff
  file 'no_hypothetical' from noHypothetical_gff
  file 'hypothetical' from hypothetical_gff
  file 'transposase' from transposase_gff
  file 'cpg' from cpg_bedGraph
  file 'gpc' from gpc_bedGraph
  file 'dam' from dam_bedGraph
  file 'dcm' from dcm_bedGraph
  file 'chr.sizes' from chr_sizes

  output:
  file "*" optional true

  """
  # Get Files
  cp -R /work/jbrowse/* . ;

  # Format FASTA file for JBROWSE
  prepare-refseqs.pl --fasta $input --key \"${params.prefix}\" --out \"data\" ;

  # Add GC content Track
  bedGraphToBigWig GC_content.bedGraph GC_content.sizes data/GC_content.bw ;
  add-bw-track.pl --bw_url GC_content.bw --plot --label "GC Content" --key "GC Content" \
  --category "GC Content" --pos_color darkgray ;

  # Add track with all features
  flatfile-to-json.pl --gff $gff --key \"All features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} annotated features\" --out \"data\" ;

  # Add tRNA track
  awk '{ if (\$3 == "tRNA" ) print }' ${gff} > tRNAs.gff ;
  [ ! -s tRNAs.gff ] || flatfile-to-json.pl --gff tRNAs.gff --key \"tRNA Sequences\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} tRNA sequences\" \
  --config '{ "style": { "color": "darkgreen" }, "displayMode": "compact" }' --out \"data\" ;

  remove-track.pl --trackLabel \"${params.prefix} tRNA sequences\" --dir data &> /tmp/error

  [ ! -s tRNAs.gff ] || echo \' { \"compress\" : 0, \
                                 \"displayMode\" : \"compact\", \
                                 \"key\" : \"tRNA Sequences\", \
                                 \"label\" : \"${params.prefix} tRNA sequences\", \
                                 \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                 \"style\" : { \"className\" : \"feature\", \"color\": \"darkgreen\" }, \
                                 \"trackType\" : \"CanvasFeatures\", \
                                 \"type\" : \"CanvasFeatures\", \
                                 \"urlTemplate\" : \"tracks/${params.prefix} tRNA sequences/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json


  # Add track without hypothetical features|proteins
  [ ! -s no_hypothetical ] || flatfile-to-json.pl --gff no_hypothetical --key \"Not hypothetical features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} not hypothetical features\" --out \"data\" ;

  # Add track with all hypothetical features|proteins
  [ ! -s hypothetical ] || flatfile-to-json.pl --gff hypothetical --key \"Only hypothetical features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} only hypothetical features\" --out \"data\" ;

  # Add track with all transposases
  [ ! -s transposase ] || flatfile-to-json.pl --gff transposase --key \"Transposases\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} only transposases\" --out \"data\" ;

  # Add track with virulence features
  [ ! -s virulence ] || flatfile-to-json.pl --gff virulence --key \"Virulence features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} virulence features\" --out \"data\" ;

  # Add track with resistance features
  [ ! -s resistance ] || flatfile-to-json.pl --gff resistance --key \"All resistance features from all sources\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} resistance features from all sources\" --out \"data\" ;

  # Add track with resistance AMRFinder features
  [ ! -s amrfinder ] || flatfile-to-json.pl --gff amrfinder --key \"AMRFinderPLus features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} resistance features from AMRFinderPlus\" --out \"data\" ;

  # Add track with resistance RGI features
  [ ! -s rgi ] || flatfile-to-json.pl --gff rgi --key \"CARD-RGI features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} resistance features from CARD-RGI\" --out \"data\" ;

  # Add track with ICEs features
  [ ! -s ices ] || flatfile-to-json.pl --gff ices --key \"ICE features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} integrative and conjugative elements\" --out \"data\" ;

  # Add track with prophage features
  [ ! -s prophage ] || flatfile-to-json.pl --gff prophage --key \"Prophage features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} prophage features\" --out \"data\" ;

  # Add track with prophage sequences
  [ ! -s prophages.bed ] || flatfile-to-json.pl --bed prophages.bed --key \"Prophage Sequences\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} prophage sequences\" \
  --config '{ "style": { "color": "blue" }, "displayMode": "compact" }' --out \"data\" ;

  remove-track.pl --trackLabel \"${params.prefix} prophage sequences\" --dir data &> /tmp/error

  [ ! -s prophages.bed ] || echo \' { \"compress\" : 0, \
                                     \"displayMode\" : \"compact\", \
                                     \"key\" : \"Prophage Sequences\", \
                                     \"label\" : \"${params.prefix} prophage sequences\", \
                                     \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                     \"style\" : { \"className\" : \"feature\", \"color\": \"blue\" }, \
                                     \"trackType\" : \"CanvasFeatures\", \
                                     \"type\" : \"CanvasFeatures\", \
                                     \"urlTemplate\" : \"tracks/${params.prefix} prophage sequences/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json

  # Add track with GIs
  [ ! -s all_GIs.bed ] || flatfile-to-json.pl --bed all_GIs.bed --key \"Genomic Islands\" \
                              --trackType CanvasFeatures --trackLabel \"${params.prefix} genomic islands\" \
                              --config '{ "style": { "color": "cyan" }, "displayMode": "compact" }' --out \"data\" ;

  # Remove track for configuration
  [ ! -s all_GIs.bed ] || remove-track.pl --trackLabel \"${params.prefix} genomic islands\" --dir data &> /tmp/error

  # Re-create
  [ ! -s all_GIs.bed ] || echo \' { \"compress\" : 0, \
                                    \"displayMode\" : \"compact\", \
                                    \"key\" : \"Genomic Islands\", \
                                    \"label\" : \"${params.prefix} genomic islands\", \
                                    \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                    \"style\" : { \"className\" : \"feature\", \"color\": \"cyan\" }, \
                                    \"trackType\" : \"CanvasFeatures\", \
                                    \"type\" : \"CanvasFeatures\", \
                                    \"urlTemplate\" : \"tracks/${params.prefix} genomic islands/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json

  # Add track with rRNA sequences
  [ ! -s rrna.gff ] || flatfile-to-json.pl --gff rrna.gff --key \"rRNA Sequences\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} rRNA sequences\" \
  --config '{ "style": { "color": "blue" }, "displayMode": "compact" }' --out \"data\" ;

  remove-track.pl --trackLabel \"${params.prefix} rRNA sequences\" --dir data &> /tmp/error

  [ ! -s rrna.gff ] || echo \' { \"compress\" : 0, \
                                 \"displayMode\" : \"compact\", \
                                 \"key\" : \"rRNA Sequences\", \
                                 \"label\" : \"${params.prefix} rRNA sequences\", \
                                 \"storeClass\" : \"JBrowse/Store/SeqFeature/NCList\", \
                                 \"style\" : { \"className\" : \"feature\", \"color\": \"blue\" }, \
                                 \"trackType\" : \"CanvasFeatures\", \
                                 \"type\" : \"CanvasFeatures\", \
                                 \"urlTemplate\" : \"tracks/${params.prefix} rRNA sequences/{refseq}/trackData.json\" } \' | add-track-json.pl  data/trackList.json

  # Add track with efflux features
  [ ! -s efflux ] || flatfile-to-json.pl --gff efflux --key \"Efflux [pumps] features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} efflux features\" --out \"data\" ;

  # Add track with conjugation features
  [ ! -s conjugation ] || flatfile-to-json.pl --gff conjugation --key \"Conjugation related features\" \
  --trackType CanvasFeatures --trackLabel \"${params.prefix} Conjugation related features\" --out \"data\" ;

  # Format bedGraphs
  ## cpg
  [ ! -s cpg ] || bedGraphToBigWig cpg chr.sizes data/cpg.bw ;
  ## gpc
  [ ! -s gpc ] || bedGraphToBigWig gpc chr.sizes data/gpc.bw ;
  ## dam
  [ ! -s dam ] || bedGraphToBigWig dam chr.sizes data/dam.bw ;
  ## dcm
  [ ! -s dcm ] || bedGraphToBigWig dcm chr.sizes data/dcm.bw ;

  # Add BigWigs
  ## cpg
  [ ! -s data/cpg.bw ] || add-bw-track.pl --bw_url cpg.bw --plot --label "CpG Methylations" --key "CpG Methylations" --category "Methylations" \
  --pos_color blue ;
  ## gpc
  [ ! -s data/gpc.bw ] || add-bw-track.pl --bw_url gpc.bw --plot --label "GpC Methylations" --key "GpC Methylations" --category "Methylations" \
  --pos_color purple ;
  ## dam
  [ ! -s data/dam.bw ] || add-bw-track.pl --bw_url dam.bw --plot --label "Dam Methylations" --key "Dam Methylations" --category "Methylations" \
  --pos_color pink ;
  ## dcm
  [ ! -s data/dcm.bw ] || add-bw-track.pl --bw_url dcm.bw --plot --label "Dcm Methylations" --key "Dcm Methylations" --category "Methylations" \
  --pos_color cyan ;
  """
}

process SQL_db {
  publishDir "${outDir}/sqlDB", mode: 'copy'
  container 'fmalmeida/bacannot:renv'

  input:
  file gff from final_gff
  file genes_nt from genes_nt_sql
  file genes_aa from genes_aa_sql

  output:
  file "${prefix}.sqlite"
  val 'finished' into finish

  script:
  """
  # Convert a FASTA to tabular (bioawk)
  bioawk -c fastx '{print \$name"\\t"\$comment"\\t"\$seq}' $genes_nt > genes_nt.tsv ;
  bioawk -c fastx '{print \$name"\\t"\$comment"\\t"\$seq}' $genes_aa > genes_aa.tsv ;

  # Create SQL db
  gff2sql.R -i $gff -o ${prefix}.sqlite -n genes_nt.tsv -a genes_aa.tsv &> /tmp/log
  """
}

// Getting config file as a file
configFile = file(workflow.configFiles[0])

process report {
  publishDir "${outDir}/report_files", mode: 'copy'
  container 'fmalmeida/bacannot:renv'

  input:
  val x from finish
  file 'final.gff' from final_gff
  file rgi_table from rgi_output
  file rgi_perfect from rgi_perfect
  file rgi_strict from rgi_strict
  file amrfinder_result from amrfinder_output
  file amrfinder_summary from amrfinderplus_table
  file 'resistance.gff' from resistance_gff
  file vfdb_blast from vfdb_blast_genes2
  file victors_blast from victors_blast_genes2
  file ice_blast from ice_blast_genes2
  //file virsorter_csv from virsorter_csv.ifEmpty('virsorter_empty')
  file phigaro_txt from phigaro_txt.ifEmpty('phigaro_empty')
  file phast_blast from phast_blast_genes2

  output:
  file '*.html'

  script:
  """
  cp /work/rscripts/*.Rmd . ;

  ## Generate Resistance Report
  Rscript -e 'rmarkdown::render("report_resistance.Rmd", params = list(\
    amrfinder = "$amrfinder_result", \
    query = "${params.prefix}", \
    rgitool = "$rgi_table", \
    rgiperfect = "$rgi_perfect", \
    rgistrict = "$rgi_strict", \
    gff_resistance = "resistance.gff", \
    ncbi_amr = "${amrfinder_summary}"))'

  ## Generate Virulence Report
  Rscript -e 'rmarkdown::render("report_virulence.Rmd" , \
  params = list( vfdb_blast = "${vfdb_blast}", \
                 blast_id = ${params.diamond_virulence_identity} , \
                 blast_cov = ${params.diamond_virulence_queryCoverage},
                 gff = "final.gff",
                 victors_blast = "${victors_blast}",
                 query = "${params.prefix}"))'

  ## Generate MGEs report
  Rscript -e 'rmarkdown::render("report_MGEs.Rmd", params = list( \
                 phigaro_dir = "../prophages/phigaro",
                 phigaro_txt = "${phigaro_txt}",
                 ice_prot_blast = "${ice_blast}",
                 query = "${params.prefix}",
                 gff = "final.gff",
                 blast_id = ${params.diamond_MGEs_identity},
                 blast_cov = ${params.diamond_MGEs_queryCoverage},
                 phast_blast = "${phast_blast}"))'
  """
}

/*
                                      Set log message
*/
log.info "========================================="
log.info " Docker-based Genome Annotation Pipeline "
log.info "========================================="
def summary = [:]
summary['Input fasta']  = params.genome
summary['Output prefix']   = params.prefix
summary['Output dir']   = "${params.outDir}"
summary['Number of threads used'] = params.threads
summary['Blast % ID - Virulence Genes'] = params.diamond_virulence_identity
summary['Blast query coverage - Virulence Genes'] = params.diamond_virulence_queryCoverage
summary['Blast % ID - ICEs and Phages'] = params.diamond_MGEs_identity
summary['Blast query coverage - ICEs and Phages'] = params.diamond_MGEs_queryCoverage
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Configuration file'] = workflow.configFiles[0]
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

// Completition message
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Execution duration: $workflow.duration"
}
