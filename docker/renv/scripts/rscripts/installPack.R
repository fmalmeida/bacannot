#!/usr/bin/Rscript

install.packages("sf", repos = "https://cloud.r-project.org/", dependencies = TRUE) ;
install.packages(c("DataCombine", "glue", "docopt", "dplyr", 
                   "stringr", "DT", "ggplot2",
                   "bookdown", "plyr", "tidyr"), repos = "https://cloud.r-project.org/", dependencies = TRUE) ;


# Install bioc packages
install.packages("BiocManager", repos = "https://cloud.r-project.org/", dependencies = TRUE) ;
BiocManager::install("Rsamtools", ask = FALSE) ;
BiocManager::install("ballgown", ask = FALSE)
