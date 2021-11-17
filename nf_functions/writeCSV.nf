def write_csv(in_list) {

  return in_list.collectFile( name:"samples.csv", sort: { it[0] }, newLine: true ) {

    /*
     * checking given/available values
     */
    

    // Check for Illumina reads
    if (it.illumina) {

        if (it.illumina.size() == 1) {          // just one read

          fwd_pair = "missing_pairFWD"
          rev_pair = "missing_pairREV"
          unpaired = it.illumina[0]

        } else if (it.illumina.size() == 2) {   // two reads are given

          fwd_pair = it.illumina[0]
          rev_pair = it.illumina[1]
          unpaired = "missing_single"

        } else if (it.illumina.size() == 3) {   // three reads are given

          fwd_pair = it.illumina[0]
          rev_pair = it.illumina[1]
          unpaired = it.illumina[2]

        }
    
    } else {                                    // no reads are given
        fwd_pair = "missing_pairFWD"
        rev_pair = "missing_pairREV"
        unpaired = "missing_single"
    }

    
    
    // Check long reads
    if (it.nanopore) {                          // nanopore is given

      lr_type = "nanopore"
      lreads  = it.nanopore

    } else if (it.pacbio) {                     // pacbio is given

      lr_type = "pacbio"
      lreads  = it.pacbio

    } else {                                    // none is given

      lr_type = "missing_lr_type"
      lreads  = "missing_lreads"

    }

    // Check fast5
    fast5 = (it.fast5) ? it.fast5  : "missing_fast5"

    // Check assembly
    assembly = (it.assembly) ? it.assembly : "missing_assembly"

    // Check resfinder
    // it uses the command line param as default and overwrites with
    // sample specific value if found inside samplesheet
    resfinder = (params.resfinder_species) ? params.resfinder_species : "missing_resfinder"
    if (it.resfinder) { resfinder = it.resfinder  }

    // Check entrypoint
    if (assembly != "missing_assembly") {

      // assembled genome is given
      "${it.id},annotation,${fwd_pair},${rev_pair},${unpaired},${lreads},${lr_type},${fast5},${assembly},${resfinder}"

    } else if (assembly == "missing_assembly" && 
               fwd_pair == "missing_pairFWD"  && 
               rev_pair == "missing_pairREV"  && 
               unpaired == "missing_single"   && 
               lreads   != "missing_lreads") {

      // short reads are not available but long reads are
      "${it.id},flye,${fwd_pair},${rev_pair},${unpaired},${lreads},${lr_type},${fast5},${assembly},${resfinder}"

    } else {

      // short reads are available
      "${it.id},unicycler,${fwd_pair},${rev_pair},${unpaired},${lreads},${lr_type},${fast5},${assembly},${resfinder}"

    }
  }

}
