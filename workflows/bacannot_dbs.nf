/*
 * Include modules
 */
include { PROKKA_DB        } from '../modules/bacannot_dbs/prokka.nf'
include { MLST_DB          } from '../modules/bacannot_dbs/mlst.nf'
include { CARD_DB          } from '../modules/bacannot_dbs/card.nf'
include { PLATON_DB        } from '../modules/bacannot_dbs/platon.nf'
include { RESFINDER_DB     } from '../modules/bacannot_dbs/resfinder.nf'
include { PLASMIDFINDER_DB } from '../modules/bacannot_dbs/plasmidfinder.nf'
include { PHIGARO_DB       } from '../modules/bacannot_dbs/phigaro.nf'
include { AMRFINDER_DB     } from '../modules/bacannot_dbs/amrfinder.nf'
include { ARGMINER_DB      } from '../modules/bacannot_dbs/argminer.nf'
include { VFDB_DB          } from '../modules/bacannot_dbs/vfdb.nf'
include { VICTORS_DB       } from '../modules/bacannot_dbs/victors.nf'
include { ICEBERG_DB       } from '../modules/bacannot_dbs/iceberg.nf'
include { PHAST_DB         } from '../modules/bacannot_dbs/phast.nf'
include { KOFAMSCAN_DB     } from '../modules/bacannot_dbs/kofamscan.nf'

/*
    DEF WORKFLOW
*/

workflow CREATE_DBS {

    // prokka database
    // core pipeline -- cannot skip
    download_db("prokka", "PROKKA_DB")

    // mlst database
    // core pipeline -- cannot skip
    download_db("mlst", "MLST_DB")

    /*
     * kofamscan -- can skip
     */
    if (!params.skip_kofamscan_db) { download_db("kofamscan", "KOFAMSCAN_DB") }

    /*
     * resistance -- can skip
     */
    // card database
    if (!params.skip_card_db) { download_db("card", "CARD_DB") }

    // resfinder database
    if (!params.skip_resfinder_db) { download_db("resfinder", "RESFINDER_DB") }

    // amrfinder database
    if (!params.skip_amrfinder_db) { download_db("amrfinder", "AMRFINDER_DB") }

    // argminer database
    if (!params.skip_argminer_db) { download_db("argminer", "ARGMINER_DB") }

    /*
     * Plasmids -- can skip
     */
    // platon database
    if (!params.skip_platon_db) { download_db("platon", "PLATON_DB") }    

    // plasmidfinder database
    if (!params.skip_plasmidfinder_db) { download_db("plasmidfinder", "PLASMIDFINDER_DB") }

    /*
     * Prophages -- can skip
     */
    // phigaro database
    if (!params.skip_phigaro_db) { download_db("phigaro", "PHIGARO_DB") }

    // phast database
    if (!params.skip_phast_db) { download_db("phast", "PHAST_DB") }

    /*
     * Virulence -- can skip
     */
    // vfdb database
    if (!params.skip_vfdb_db) { download_db("vfdb", "VFDB_DB") }

    // victors database
    if (!params.skip_victors_db) { download_db("victors", "VICTORS_DB") }

    /*
     * ICEs -- can skip
     */
    // iceberg database
    if (!params.skip_iceberg_db) { download_db("iceberg", "ICEBERG_DB") }

}

def download_db(database, module) {
    if (file("${params.output}/${database}_db").exists() && params.force_update == false) {
        println "NOTE:\n\t=> ${database} database already exists and --force_update was not used. Skipping."
    } else {
        "${module}"()
    } 
}