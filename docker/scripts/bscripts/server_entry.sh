#!/bin/bash
sequenceserver &> /tmp/sequenceserver.log &
Rscript -e 'rmarkdown::run("shiny_parser.Rmd", shiny_args = list(launch.browser = T, host = "0.0.0.0", port = 3838))'