#!/usr/bin/Rscript
# Setting help
'usage: correctSourceAndFeature.R [--input=<file> --out=<chr>]

options:
  -i, --input=<file>    GFF to add annotated features source and type based on Pfam subset.
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

# Set function
# This function is used to extract the values of the fields stored in
# Attributes column of gff file.
getmotif <- function (x, field, attrsep = ";") { 
  s = strsplit(x, split = attrsep, fixed = TRUE) 
  sapply(s, function(atts) { 
    a = strsplit(atts, split = "=", fixed = TRUE) 
    m = match(field, sapply(a, "[", 1)) 
    if (!is.na(m)) { rv = a[[m]][2] 
    } 
    else { 
      rv = as.character(NA) 
    } 
    return(rv) 
  }) 
}

# Function used to remove redundancy
reduce_row = function(i) {
  d <- unlist(strsplit(i, split=","))
  paste(unique(d), collapse = ',') 
}

# Read gff file
gff <- gffRead(opt$input)
gff$attributes <- gsub(" ", "_", gff$attributes)
gff$attributes <- gsub("protein motif:", "protein_motif=", gff$attributes)
gff$attributes <- gsub(",protein_motif:", ";protein_motif=", gff$attributes)

# First step is to subset those entries that have something related to one 
# of the pfam subsets

sub <- grepl.sub(gff, pattern = "_subset", Var = "attributes")
not <- grepl.sub(gff, pattern = "_subset", Var = "attributes", keep.found = FALSE)

# Secondly we need to store the previous value of the source column in 
# order to add to it the name of the pfam subset database

previous_source <- sub$source

# Then, we need to extract the name of the database 
# and feature type.
motifs <- getmotif(sub$attributes, "protein_motif", ";")
sources <- sapply(strsplit(motifs, split = "_", fixed = TRUE), "[", 1)
features <- sapply(strsplit(motifs, split = "_", fixed = TRUE), "[", 2)

# Finally, we add these names to the previous sources
new_source <- paste(previous_source, sources, sep = ",")
sub$source <- new_source

# Then, we add it to previously written features
previous_features <- sub$feature
new_features <- paste(previous_features, features, sep = ",")
sub$feature <- new_features

# Merge gff subsets
merged_df <- merge.data.frame(sub, not, all = TRUE)
feat <- merged_df$feature
merged_df$feature <- sapply(feat, reduce_row)
source <- merged_df$source
merged_df$source <- sapply(source, reduce_row)
merged_df <- merged_df[order(merged_df$seqname, merged_df$start),]

# Write output
write.table(merged_df, file = opt$out, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE)
