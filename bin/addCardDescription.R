#!/usr/bin/Rscript
# Setting help
'usage: addCardDescription.R [--input=<file> --gff=<file> --out=<chr> --database=<chr> --type=<chr> --scoverage=<int>]

options:
  -i, --input=<file>    Tabular Blast to be added to GFF
  -g, --gff=<file>      GFF file to add Blast hits into
  -o, --out=<chr>       Output file name [default: out.gff]
  -d, --database=<chr>  Name of databased which Blast came from
  -t, --type=<chr>      Type of feature blasted. Ex: resistance
  -c, --scoverage=<int> Minimum subject coverage to keep' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)
if (is.null(opt$input)){
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

# Load libraries
suppressMessages(library(DataCombine))
suppressMessages(library(ballgown))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

# Load CARD entries index. These will be used to write
# the attributes columns of CARD entries.
cat_index <- read.table("/work/indexes/aro_categories_index.csv", header = TRUE, sep = "\t")
cat <- read.table("/work/indexes/aro_categories.csv", header = TRUE, sep = "\t")
index <- read.table("/work/indexes/aro_index.csv", header = TRUE, sep = "\t", fill = TRUE)

# Function used to remove redundancy
reduce_row = function(i) {
  d <- unlist(strsplit(i, split=","))
  paste(unique(d), collapse = ',') 
}

getAttributeField <- function (x, field, attrsep = ";") { 
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

# Check if file is empty
if (file.info(opt$input)$size > 0 ) {
# Merge indexes to create a full index with all CARD values
merged <- merge.data.frame(cat_index, index, by.x = "Protein.Accession",
                           by.y = "Protein.Accession", all = TRUE)
card_indexes <- merge.data.frame(merged, cat, by.x = "ARO.Accession", 
                                 by.y = "ARO.Accession", all = TRUE)

# Load Blast tabular file
blastFile <- read.delim(opt$input, header = FALSE)
blastHeader <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
                 "qend", "sstart", "send", "slen", "evalue", "bitscore", "stitle")
colnames(blastFile) <- blastHeader

# Filter blast based on subject coverage
if (!is.null(opt$scoverage)) {
blastFile$scov <- (blastFile$length / blastFile$slen) * 100
blastFile <- dplyr::filter(blastFile, scov >= as.integer(opt$scoverage))}

if (nrow(blastFile) > 0) {
# Remove duplicates based on bitscore
blastFile <- blastFile[order(blastFile$qseqid, -abs(blastFile$pident), -abs(blastFile$scov), -abs(blastFile$bitscore) ), ]
blastFile <-blastFile[ !duplicated(blastFile$qseqid), ]
blastFile <- blastFile[order(blastFile$qseqid),]

## Filter per %Identity
# blast_filtered <- subset(blastFile, pident >= as.integer(opt$pident))
blast_filtered <- blastFile
ssids <- as.vector(blast_filtered$sseqid)
aroID <- sapply(strsplit(ssids, "\\|"), `[`, 3)
blast_filtered$ARO <- aroID

# Load gff file to merge entries
gff <- gffRead(opt$gff)

# Create a column in gff with ids
gff$ID <- getAttributeField(gff$attributes, "ID", ";")

## Subset Card indexes
card_subset <- grepl.sub(card_indexes, pattern = aroID, Var = "ARO.Accession")
card_subset$Drug.Class <- gsub(";", ":", card_subset$Drug.Class)
card_subset$Drug.Class <- gsub(" ", "_", card_subset$Drug.Class)
card_subset$Drug.Class <- gsub("\t", "_", card_subset$Drug.Class)
card_subset$Drug.Class <- gsub("-", "_", card_subset$Drug.Class)
card_subset$AMR.Gene.Family <- gsub(";", ":", card_subset$AMR.Gene.Family)
card_subset$AMR.Gene.Family <- gsub(" ", "_", card_subset$AMR.Gene.Family)
card_subset$AMR.Gene.Family <- gsub("\t", "_", card_subset$AMR.Gene.Family)
card_subset$AMR.Gene.Family <- gsub("-", "_", card_subset$AMR.Gene.Family)
card_subset$Resistance.Mechanism <- gsub(";", ":", card_subset$Resistance.Mechanism)
card_subset$Resistance.Mechanism <- gsub(" ", "_", card_subset$Resistance.Mechanism)
card_subset$Resistance.Mechanism <- gsub("-", "_", card_subset$Resistance.Mechanism)
card_subset$Resistance.Mechanism <- gsub("\t", "_", card_subset$Resistance.Mechanism)
card_subset$Model.Name <- gsub(" ", "_", card_subset$Model.Name)
card_subset$Model.Name <- gsub("\t", "_", card_subset$Model.Name)
card_subset$Model.Name <- gsub("-", "_", card_subset$Model.Name)

# Get desired values for attributes columns
description <- paste("Additional_database=", opt$database, ";", 
                     "ARO=", card_subset$ARO.Accession, ";", "Gene_Family=", 
                     card_subset$AMR.Gene.Family, ";", "Drug_Class=", card_subset$Drug.Class, 
                     ";", "Resistance_Mechanism=", card_subset$Resistance.Mechanism, ";", 
                     "DB_Name=", card_subset$Model.Name, ";CVTERM=", card_subset$CVTERM, sep = "")

card_subset$CARD_attributes <- description
card_subset$CARD_attributes <- gsub(" ", "_", card_subset$CARD_attributes)
card_subset$CARD_attributes <- gsub("\t", "_", card_subset$CARD_attributes)

# Concatenate new attributes values
blast_filtered <- merge.data.frame(blast_filtered, card_subset, by.x = "ARO", 
                                   by.y = "ARO.Accession", all = TRUE)

blast_filtered <- blast_filtered[order(blastFile$qseqid),]

# Get gene names from blast hits
ids <- blast_filtered$qseqid

# Subset based on gene IDs
sub <- grepl.sub(gff, pattern = ids, Var = "ID") %>% select(seqname, source, feature, start, end, score, strand, frame, attributes, ID)
not <- grepl.sub(gff, pattern = ids, Var = "ID", keep.found = FALSE) %>% select(seqname, source, feature, start, end, score, strand, frame, attributes)

# Change fields - Add database source and feature type
## source
s <- sub$source
sn <- opt$database
snew <- paste(s, sn, sep = ",")
sub$source <- snew

## feature
f <- sub$feature
fn <- opt$type
fnew <- paste(f, fn, sep = ",")
sub$feature <- fnew

## attributes
sub <- merge.data.frame(sub, blast_filtered, by.x = "ID", 
                        by.y = "qseqid", all = TRUE)
sub <- unite(sub, "attributes", c("attributes", "CARD_attributes"), sep = ";") %>%
  select(seqname, source, feature, start, end, score, strand, frame, attributes)

# Merge files
merged_df <- merge.data.frame(sub, not, all = TRUE)
feat <- merged_df$feature
merged_df$feature <- sapply(feat, reduce_row)
source <- merged_df$source
merged_df$source <- sapply(source, reduce_row)
merged_df <- merged_df[order(merged_df$seqname, merged_df$start),]

# Write output
write.table(merged_df, file = opt$out, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE, append = FALSE)
} else {
  # Load GFF file
  gff <- gffRead(opt$gff)
  # Write output
  write.table(gff, file = opt$out, quote = FALSE, sep = "\t", 
              col.names = FALSE, row.names = FALSE)
}} else {
  # Load GFF file
  gff <- gffRead(opt$gff)
  # Write output
  write.table(gff, file = opt$out, quote = FALSE, sep = "\t", 
              col.names = FALSE, row.names = FALSE)
}
