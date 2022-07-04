//
// This file holds several functions specific to the workflow/bacannot.nf in the fmalmeida/bacannot pipeline
//

class WorkflowBacannot {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        // input has been given and user does not want to download databases?
        if (!params.input && !params.get_dbs) {
            log.error "Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.yml'. Or select the download databases mode with --get_dbs."
            System.exit(1)
        }

        // prokka params checkup
        if (params.prokka_kingdom && !params.prokka_genetic_code) {
            log.error """
            ERROR!

            A minor error has occurred
            ==> User have set --prokka_kingdom but forgot --prokka_genetic_code.

            These parameters must be used together. If you change prokka defaults kingdom parameter you must set the genetic code to be used for translation.

            If in doubt with these parameters let it blank, or get more information in Prokka's documentation.

            Cheers.
            """.stripIndent()
            System.exit(1)
        }

    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

}
