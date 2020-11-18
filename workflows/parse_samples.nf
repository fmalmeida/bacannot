workflow parse_samplesheet {

  take:
    data
  main:

    // Parse input to check for missing entries
    parsed = []

    data.each {

      // Check for Illumina reads
      if (it.illumina) {

        // Only one read
        if (it.illumina.size() == 1) {
          it['paired'] = "missing_paired"
          it['single'] = it.illumina[0]
        } else if (it.illumina.size() == 2) {
          it['paired'] = [it.illumina[0], it.illumina[1]]
          it['single'] = "missing_single"
        } else if (it.illumina.size() == 3) {
          it['paired'] = [it.illumina[0], it.illumina[1]]
          it['single'] = it.illumina[2]
        }
      } else {
        it['paired'] = "missing_paired"
        it['single'] = "missing_single"
      }

      // Check long reads
      if (it.nanopore) {
        it['lr_type'] = "nanopore"
        it['lreads'] = it.nanopore
      } else if (it.pacbio) {
        it['lr_type'] = "pacbio"
        it['lreads'] = it.pacbio
      } else {
        it['lr_type'] = "missing_lr_type"
        it['lreads'] = "missing_lreads"
      }

      // Check fast5
      it['fast5'] = (it.fast5)   ? it.fast5  : "missing_fast5"

      // Check assembly
      it['assembly'] = (it.assembly) ? it.assembly : "missing_assembly"

      // Save
      parsed.add(it)
    }

    emit:
      Channel.from(parsed)

}
