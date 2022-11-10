#!/usr/bin/env Rscript
# A supporting script for writing amrfinderplus results as csv!
# To be used with fmalmeida/bacannot-compare
# Author: Felipe Marques de Almeida <almeidafmarques@gmail.com>

#################################################
### Setting Help message and input parameters ###
#################################################
'usage: amrfinder2csv.R [ --input=<file> --output=<chr> --aro=<file> --sample=<chr> ]
options:
-i, --input=<file>       Bacannot JSON summary
-o, --output=<chr>       Output file [default: out.tsv]
-a, --aro=<file>         CARD ARO categories
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
# functions for CARD metadata
#
aro_index <- 
    read.csv( opt$aro, sep = "\t" )
aro_index$ARO.Accession <- 
    str_remove(aro_index$ARO.Accession, "ARO:")

get_aro_category <- function(aro) {
    #aro <- 3003923
    df <- aro_index %>% dplyr::filter(ARO.Accession == aro)
    result <- str_split(df$Drug.Class, ";")[[1]][1] # tentar padronizar pegando so o primeiro
    return(result)
}

get_aro_gene_name <- function(aro) {
    # aro <- 3003923
    df <- aro_index %>% filter(ARO.Accession == aro)
    result <- df$CARD.Short.Name
    return(result)
}

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
    without_aro      <- 
    setNames(
        data.frame(
        matrix(ncol = 6, nrow = 0)
        ), c("sample", "contig", "start", "end", "gene", "subclass")
    )

#
# parse AMRFinderPlus Results
#
for (sample in samples) {
    subset <- summary_file[[sample]]$resistance$amrfinderplus
    all_entries <- names(subset)
    for (entry in all_entries) {
        if (entry != "total") {
            
            if (subset[[entry]][['type']] == 'AMR') {
                
                contig = subset[[entry]]$contig
                start  = subset[[entry]]$start
                end    = subset[[entry]]$end
                
                if (is.null(subset[[entry]][['card_aro']])) {
                subclass = subset[[entry]]$subclass
                without_aro[nrow(without_aro) + 1,] <- 
                    c( sample, contig, start, end, subset[[entry]]$gene, subclass )
                } else {
                if (is.null(subset[[entry]][['subclass']])) {
                    subclass = "None"
                } else {
                    subclass = get_aro_category(subset[[entry]]$card_aro)
                }
                identified_genes[nrow(identified_genes) + 1,] <- 
                    c( sample, contig, start, end, subset[[entry]]$gene, subclass )  
                }
                
            }
        }
    } 
}

#
# filter results
#
if (as.character(opt$sample) != 'all') {
    identified_genes <-
        identified_genes %>% dplyr::filter( as.character(sample) == as.character(opt$sample) )
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