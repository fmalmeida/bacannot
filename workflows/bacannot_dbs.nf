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

}

def download_db(database, module) {
    if (file("${params.output}/${database}_db").exists() && params.force_update == false) {
        println "NOTE:\n\t=> ${database} database already exists and --force_update was not used. Skipping."
    } else {
        "${module}"()
    } 
}