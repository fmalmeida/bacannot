#!/usr/bin/Rscript

# Setting parameters
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}


df <- read.delim(opt$input)
df$Contig <- sub("_[^_]+$", "", df$Contig)
df$source <- "RGI"
df$type <- "resistance"
df$frame <- "."
df$score <- "."

attributes <- paste0("ID=", df$ID, ";", "gene=", df$Best_Hit_ARO, ";", "product=", df$Drug.Class, ";", "gene_family=",
                     df$AMR.Gene.Family, ";", "resistance_mechanism=", df$Resistance.Mechanism)
df$attributes <- attributes

gff <- df[, c("Contig", "source", "type", "Start", "Stop", "score", "Orientation", "frame", "attributes")]

colnames(gff) <- c("seqname", "source", "type", "start", "end", "score", "strand", "frame", "attributes")


##Write file

write.table(gff, file = opt$out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
