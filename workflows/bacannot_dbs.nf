/*
 * Include modules
 */
include { CARD_DB          } from '../modules/bacannot_dbs/card.nf'
include { PLATON_DB        } from '../modules/bacannot_dbs/platon.nf'
include { RESFINDER_DB     } from '../modules/bacannot_dbs/resfinder.nf'
include { PLASMIDFINDER_DB } from '../modules/bacannot_dbs/plasmidfinder.nf'

/*
    DEF WORKFLOW
*/

workflow CREATE_DBS {

    // card database
    if (!params.skip_card_db) { download_db("card", "CARD_DB") }

    // platon database
    if (!params.skip_platon_db) { download_db("platon", "PLATON_DB") }

    // resfinder database
    if (!params.skip_resfinder_db) { download_db("resfinder", "RESFINDER_DB") }

    // plasmidfinder database
    if (!params.skip_plasmidfinder_db) { download_db("plasmidfinder", "PLASMIDFINDER_DB") }

}

def download_db(database, module) {
    if (file("${params.output}/${database}_db").exists() && params.force_update == false) {
        println "NOTE:\n\t=> ${database} database already exists and --force_update was not used. Skipping."
    } else {
        "${module}"()
    } 
}