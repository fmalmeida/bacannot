/*
 * Include modules (Execution setup)
 */

include { MAKE_KARYOTYPE } from '../modules/generic/karyotype'
include { GC_SKEW        } from '../modules/generic/gc_skew'
include { AMRFINDER2TSV  } from '../modules/resistance/amrfinder2tsv'

/*
    DEF WORKFLOW
*/

workflow CIRCOS {
    take:
        input_ch
        conf_ch
    
    main:

    //
    // generate karyotype files
    //
    MAKE_KARYOTYPE( input_ch ).karyotype.view()

    //
    // calculate gc skew
    //
    GC_SKEW( input_ch ).skew.view()

    //
    // convert amrfinder annotation to tsv
    //
    AMRFINDER2TSV( 
        input_ch,
        file( "$baseDir/assets/aro_index.tsv" )
    )

}
