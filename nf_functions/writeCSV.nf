def write_csv(in_list) {

  return in_list.collectFile( name:"samples.txt", sort: { it[0] }, newLine: true ) {

    // Check short read pairs
    fwd_pair = (it.paired != "missing_paired") ? "${it.paired[0]}" : "missing_pairFWD"
    rev_pair = (it.paired != "missing_paired") ? "${it.paired[1]}" : "missing_pairREV"

    // Check lreads technology
    lr_type = (it.lr_type != "missing_lr_type") ? "${it.lr_type}" : "missing_lr_type"

    // Check entrypoint
    if (it.assembly != "missing_assembly") {

      "${it.id},annotation,${fwd_pair},${rev_pair},${it.single},${it.lreads},${lr_type},${it.fast5},${it.assembly},${it.resfinder}"

    } else if (it.assembly == "missing_assembly" && it.paired == "missing_paired" && it.single == "missing_single" && it.lreads != "missing_lreads") {

      "${it.id},flye,${fwd_pair},${rev_pair},${it.single},${it.lreads},${lr_type},${it.fast5},${it.assembly},${it.resfinder}"

    } else {

      "${it.id},unicycler,${fwd_pair},${rev_pair},${it.single},${it.lreads},${lr_type},${it.fast5},${it.assembly},${it.resfinder}"

    }
  }

}
