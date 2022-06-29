include { write_csv } from '../nf_functions/writeCSV.nf'
workflow PARSE_SAMPLESHEET {

  take:
    data
  main:

    // iterate over input list
    custom_csv_ch = write_csv(Channel.fromList(data))
  
    // now we parse the csv created
    parsed_csv_ch = custom_csv_ch.splitCsv(header: ['name', 'entrypoint', 'fwd', 'rev', 'single', 'lreads', 'lr_type', 'fast5', 'assembly', 'resfinder']).map{ row ->
      tuple(
        row.name, 
        row.entrypoint,
        (row.fwd == "missing_pairFWD") ? row.fwd : file(row.fwd, checkIfExists: true), 
        (row.rev == "missing_pairREV") ? row.rev : file(row.rev, checkIfExists: true),
        (row.single == "missing_single") ? row.single : file(row.single, checkIfExists: true), 
        (row.lreads == "missing_lreads") ? row.lreads : file(row.lreads, checkIfExists: true), 
        row.lr_type,
        (row.fast5 == "missing_fast5") ? row.fast5 : file(row.fast5, checkIfExists: true), 
        (row.assembly == "missing_assembly") ? row.assembly : file(row.assembly, checkIfExists: true), 
        row.resfinder
      )
    }

    emit:
    parsed_csv_ch

}
