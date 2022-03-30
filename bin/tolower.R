#!/usr/bin/Rscript
# Setting Help
'usage: tolower.R [--input=<file> --out=<chr> ]

options:
  -i, --input=<file>    Tabular Blast to be added to GFF
  -o, --out=<chr>       Output file name [default: out.gff]' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)

if (is.null(opt$input)){
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

# Load libraries
suppressMessages(library(ballgown))
suppressMessages(library(DataCombine))
suppressMessages(library(dplyr))

# Load GFF file
gff <- gffRead(opt$input)
# Lower case the attributes column
gff$attributes <- sapply(gff$attributes, tolower)
# Write output
write.table(gff, file = opt$out, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE)