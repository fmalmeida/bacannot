#!/usr/bin/env Rscript
# A supporting script for writing VFDB results as csv!
# To be used with fmalmeida/bacannot-compare
# Author: Felipe Marques de Almeida <almeidafmarques@gmail.com>

#################################################
### Setting Help message and input parameters ###
#################################################
'usage: vfdb2csv.R [ --input=<file> --sample=<chr> --output=<chr> ]
options:
-i, --input=<file>       Bacannot JSON summary
-o, --output=<chr>       Output file [default: out.tsv]
-s, --sample=<chr>       Filter results of only this sample [default: all]
' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)

if (is.null(opt$input)){
    stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

#
# load libraries
#
library(jsonlite)
library(stringr)
library(dplyr)

#
# read annotation summary file
#
summary_file <-
    read_json( opt$input )
samples <- 
    names( summary_file ) # get available samples in file

#
# start empty data.frames for parsing
#
identified_genes <- 
    setNames(
        data.frame(
        matrix(ncol = 6, nrow = 0)
        ), c("sample", "contig", "start", "end", "gene", "subclass")
    )

#
# parse VFDB Results
#
for (sample in samples) {
    subset <- summary_file[[sample]]$virulence$VFDB
    all_entries <- names(subset)
    for (entry in all_entries) {
        if (entry != "total") {
        
        contig   = subset[[entry]]$chr
        start    = subset[[entry]]$start
        end      = subset[[entry]]$end
        subclass = subset[[entry]]$virulence_factor
        gene     = subset[[entry]]$gene

        identified_genes[nrow(identified_genes) + 1,] <- 
        c( sample, contig, start, end, gene, subclass )
        
        }
    }
}

identified_genes <- 
    identified_genes %>% 
    distinct()

#
# filter results
#
if (as.character(opt$sample) != 'all') {
    identified_genes <-
        identified_genes %>% 
        dplyr::filter( as.character(sample) == as.character(opt$sample) )
}

#
# save table as file
#
write.table(
    identified_genes,
    file = as.character( opt$output ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
)