def parse_csv(in_ch) {

  return in_ch.splitCsv(header: ['name', 'entrypoint', 'fwd', 'rev', 'single', 'lreads', 'lr_type', 'fast5', 'assembly', 'resfinder']).map{ row ->

    tuple(
      row.name, 
      row.entrypoint,
      (row.fwd == "missing_pairFWD") ? row.fwd : file(row.fwd), 
      (row.rev == "missing_pairREV") ? row.rev : file(row.rev),
      (row.single == "missing_single") ? row.single : file(row.single), 
      (row.lreads == "missing_lreads") ? row.lreads : file(row.lreads), 
      row.lr_type,
      (row.fast5 == "missing_fast5") ? row.fast5 : file(row.fast5), 
      (row.assembly == "missing_assembly") ? row.assembly : file(row.assembly), 
      row.resfinder
    )
  
  }

}

def filter_ch(in_ch, entrypoint) {
  parse_csv(in_ch.map { it.text }) | filter { it[1] == "${entrypoint}" }
}
