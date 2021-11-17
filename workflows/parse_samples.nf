include { write_csv } from '../nf_functions/writeCSV.nf'
workflow parse_samplesheet {

  take:
    data
  main:

    // iterate over input list
    custom_csv = write_csv(Channel.fromList(data))
  
    // now we parse the csv created
    parsed_csv = custom_csv.splitCsv(header: ['name', 'entrypoint', 'fwd', 'rev', 'single', 'lreads', 'lr_type', 'fast5', 'assembly', 'resfinder']).map{ row ->
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

    emit:
    parsed_csv

}
