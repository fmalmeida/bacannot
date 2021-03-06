def parse_csv(in_ch) {

  return in_ch.splitCsv(header: ['name', 'entrypoint', 'fwd', 'rev', 'single', 'lreads', 'lr_type', 'fast5', 'assembly', 'resfinder']).map{ row ->

    if (row.entrypoint == "flye" && row.lreads != "missing_lreads" && row.lr_type == "missing_lr_type") {
      println """
      ERROR!

      A minor error has occurred
        ==> User used --lreads but forgot --lreads_type.

      When giving longreads as input, you must tell the pipeline from wich tech it comes from: 'nanopore' or 'pacbio'

      Cheers.
      """.stripIndent()

      exit 1
    } else if (row.resfinder == "other" || row.resfinder == "Other" || row.resfinder == "OTHER") {
      println """
      ERROR!

      A minor error has occurred
        ==> User has set the resfinder panel to "Other" (Genome: ${row.name})

      This is impossible, since the pipeline tries to annotation point finder mutation and these are incompatible with the "Other" panel

      Cheers.
      """.stripIndent()

      exit 1
    } else {
      tuple(row.name, row.entrypoint, (row.fwd == "missing_pairFWD") ? row.fwd : file(row.fwd), (row.rev == "missing_pairREV") ? row.rev : file(row.rev),
      (row.single == "missing_single") ? row.single : file(row.single), (row.lreads == "missing_lreads") ? row.lreads : file(row.lreads), row.lr_type,
      (row.fast5 == "missing_fast5") ? row.fast5 : file(row.fast5), (row.assembly == "missing_assembly") ? row.assembly : file(row.assembly), row.resfinder)
    }
  }

}

def filter_ch(in_ch, entrypoint) {
  parse_csv(in_ch.map { it.text }) | filter { it[1] == "${entrypoint}" }
}
