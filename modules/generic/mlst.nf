process MLST {
   publishDir "${params.output}/${prefix}", mode: 'copy', saveAs: { filename ->
     if (filename.indexOf("_version.txt") > 0) "tools_versioning/$filename"
     else "MLST/$filename"
   }
   tag "${prefix}"

   input:
   tuple val(prefix), file(genome)
   file(bacannot_db)

   output:
   tuple val(prefix), path("${prefix}_mlst_analysis.txt")  , emit: mlst optional true
   tuple val(prefix), path("${prefix}_novel_alleles.fasta"), emit: novelAlleles optional true
   path('mlst_version.txt'), emit: version

   script:
   """
   # update tool database
   mlst_dir=\$(which mlst | sed 's/bin\\/mlst//g')
   cp ${bacannot_db}/mlst_db/* -r \${mlst_dir}/db/pubmlst/
   ( cd \$mlst_dir/scripts && ./mlst-make_blast_db )

   # Save mlst tool version
   mlst --version > mlst_version.txt ;

   # run mlst
   mlst \\
       --quiet \\
       --novel ${prefix}_novel_alleles.fasta \\
       $genome > ${prefix}_mlst_analysis.txt
   """
}
