#!/bin/bash
yes '' | sequenceserver -d /work/SequenceServerDBs -H 0.0.0.0 &> /tmp/sequenceserver.log &
Rscript -e 'rmarkdown::run("shiny_parser.Rmd", shiny_args = list(launch.browser = T, host = "0.0.0.0", port = 3838))'