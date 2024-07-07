/*
 * Include modules
 */
include { PROKKA_DB        } from '../modules/bacannot_dbs/prokka'
include { MLST_DB          } from '../modules/bacannot_dbs/mlst'
include { CARD_DB          } from '../modules/bacannot_dbs/card'
include { PLATON_DB        } from '../modules/bacannot_dbs/platon'
include { RESFINDER_DB     } from '../modules/bacannot_dbs/resfinder'
include { PLASMIDFINDER_DB } from '../modules/bacannot_dbs/plasmidfinder'
include { PHIGARO_DB       } from '../modules/bacannot_dbs/phigaro'
include { AMRFINDER_DB     } from '../modules/bacannot_dbs/amrfinder'
include { ARGMINER_DB      } from '../modules/bacannot_dbs/argminer'
include { VFDB_DB          } from '../modules/bacannot_dbs/vfdb'
include { VICTORS_DB       } from '../modules/bacannot_dbs/victors'
include { ICEBERG_DB       } from '../modules/bacannot_dbs/iceberg'
include { PHAST_DB         } from '../modules/bacannot_dbs/phast'
include { KOFAMSCAN_DB     } from '../modules/bacannot_dbs/kofamscan'
include { ANTISMASH_DB     } from '../modules/bacannot_dbs/antismash'
include { GET_ZENODO_DB    } from '../modules/bacannot_dbs/get_zenodo'
include { SOURMASH_DB      } from '../modules/bacannot_dbs/sourmash'

/*
    DEF WORKFLOW
*/

workflow CREATE_DBS {

    if ( params.get_dbs && !params.get_zenodo_db ) {
        download_db("prokka", "PROKKA_DB")
        download_db("mlst", "MLST_DB")
        download_db("kofamscan", "KOFAMSCAN_DB")
        download_db("card", "CARD_DB")
        download_db("resfinder", "RESFINDER_DB")
        download_db("amrfinder", "AMRFINDER_DB")
        download_db("argminer", "ARGMINER_DB")
        download_db("platon", "PLATON_DB")
        download_db("plasmidfinder", "PLASMIDFINDER_DB")
        download_db("phigaro", "PHIGARO_DB")
        download_db("phast", "PHAST_DB")
        download_db("vfdb", "VFDB_DB")
        download_db("victors", "VICTORS_DB")
        download_db("iceberg", "ICEBERG_DB")
        download_db("antismash", "ANTISMASH_DB")
        download_db("sourmash", "SOURMASH_DB")
    } else if ( !params.get_dbs && params.get_zenodo_db ) {
        GET_ZENODO_DB()
    }

}

def download_db(database, module) {
    if (file("${params.output}/${database}_db").exists() && params.force_update == false) {
        println "NOTE:\n\t=> ${database} database already exists and --force_update was not used. Skipping."
    } else {
        "${module}"()
    } 
}