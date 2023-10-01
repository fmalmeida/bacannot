//
// This file holds several functions specific to the main.nf workflow in the fmalmeida/bacannot pipeline
//

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The pipeline\n" +
            "  https://doi.org/10.12688/f1000research.139488.1\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Print help to screen if required
    //
    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} --input bacannot_samplesheet.yaml --output ./results"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        help_string += '\n' + citation(workflow) + '\n'
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return help_string
    }

    //
    // Print parameter summary log to screen
    //
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return summary_log
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Download template config
        if (params.get_config) {
            new File("bacannot.config").write(new URL ("https://github.com/fmalmeida/bacannot/raw/master/conf/defaults.config").getText())
            log.info """
            bacannot.config file saved in working directory
            After configuration, run:
              nextflow run fmalmeida/bacannot -c ./bacannot.config
            Nice code
            """.stripIndent()
            System.exit(0)
        }

        // Download template samplesheet
        if (params.get_samplesheet) {
            new File("bacannot_samplesheet.yaml").write(new URL ("https://github.com/fmalmeida/bacannot/raw/master/example_samplesheet.yaml").getText())
            log.info """
            Samplesheet (bacannot_samplesheet.yml) file saved in working directory
            Nice code!
            """.stripIndent()
            System.exit(0)
        }

        // Download docker config
        if (params.get_docker_config) {
            new File("docker.config").write(new URL ("https://github.com/fmalmeida/bacannot/raw/master/conf/docker.config").getText())
            log.info """
            docker.config file saved in working directory
            After configuration, run:
              nextflow run fmalmeida/bacannot -c ./docker.config
            Nice code
            """.stripIndent()
            System.exit(0)
        }

        // Download singularity config
        if (params.get_singularity_config) {
            new File("singularity.config").write(new URL ("https://github.com/fmalmeida/bacannot/raw/master/conf/singularity.config").getText())
            log.info """
            singularity.config file saved in working directory
            After configuration, run:
              nextflow run fmalmeida/bacannot -c ./singularity.config
            Nice code
            """.stripIndent()
            System.exit(0)
        }

        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            NfcoreSchema.validateParameters(workflow, params, log)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        // if (params.enable_conda) {
        //     Utils.checkCondaChannels(log)
        // }

        // Check AWS batch settings
        // NfcoreTemplate.awsBatch(workflow, params)
        
    }

}
