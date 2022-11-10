/*
 * Include modules (Execution setup)
 */

include { MAKE_KARYOTYPE } from '../modules/generic/karyotype'
include { GC_SKEW        } from '../modules/generic/gc_skew'
include { AMRFINDER2TSV  } from '../modules/resistance/amrfinder2tsv'
include { VFDB2TSV       } from '../modules/virulence/vfdb2tsv'
include { PREPARE_CIRCOS } from '../modules/generic/prepare_circos'
include { PLOT           } from '../modules/generic/circos'

/*
    DEF WORKFLOW
*/

workflow CIRCOS {
    take:
        input_ch
        conf_ch
        plasmidfinder_ch
    
    main:

    //
    // generate karyotype files
    //
    MAKE_KARYOTYPE( input_ch )

    //
    // calculate gc skew
    //
    GC_SKEW( input_ch )

    //
    // convert amrfinder annotation to tsv
    //
    AMRFINDER2TSV( 
        input_ch,
        file( "$baseDir/assets/aro_index.tsv" )
    )

    //
    // convert vfdb annotation to tsv
    //
    VFDB2TSV( input_ch )

    //
    // prepare circos data
    //
    MAKE_KARYOTYPE.out.karyotype
    .join( GC_SKEW.out.skew      , remainder: true )
    .join( AMRFINDER2TSV.out.args, remainder: true )
    .join( VFDB2TSV.out.genes    , remainder: true )
    .join( plasmidfinder_ch      , remainder: true )
    .map{ it ->
        sample = it[0]
        it.remove(0)
        [ sample, it ]
    }
    .set { circos_inputs_ch }

    PREPARE_CIRCOS( circos_inputs_ch )

    //
    // plot circos
    //
    PLOT( 
        MAKE_KARYOTYPE.out.karyotype
        .join( GC_SKEW.out.skew       , remainder: true )
        .join( PREPARE_CIRCOS.out.data, remainder: true )
        .map{
            sample = it[0]
            it.remove(0)
            [ sample, it.flatten().collect() ]
        },
        conf_ch.collect()
    )

}
