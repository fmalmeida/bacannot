#!/usr/bin/Rscript
# Setting Help
'usage: addNCBIamr2Gff.R [--input=<file> --gff=<file> --out=<chr> --database=<chr> --type=<chr>]

options:
-g, --gff=<file>      GFF file to add NCBI AMR Annotations into
-i, --input=<file>    AMRFinder output
-o, --out=<chr>       Output file name [default: out.gff]
-t, --type=<chr>      Type of feature. Ex: resistance
-d, --database=<chr>  Name of databased which Blast came from' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)

if (is.null(opt$gff)){
  stop("At least one argument must be supplied (gff file)\n", call.=FALSE)
}

if (is.null(opt$input)){
  stop("At least one argument must be supplied (AMRFinder output file)\n", call.=FALSE)
}

# Load libraries
suppressMessages(library(ballgown))
suppressMessages(library(DataCombine))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))

# Function used to remove redundancy
reduce_row = function(i) {
  d <- unlist(strsplit(i, split=","))
  paste(unique(d), collapse = ',') 
}

# Function to get Attribute Fields
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
# Operator to discard patterns found
'%ni%' <- Negate('%in%')

# Load GFF File
gff <- gffRead(opt$gff)
gff$ID <- getAttributeField(as.character(gff$attributes), "ID", ";")

# Load NCBI AMRFinder output
NCBIamr <- read.delim(opt$input)

if (is.null(NCBIamr) == FALSE & dim(NCBIamr)[1] != 0) {
  
# Get its ids
ids <- NCBIamr$Protein.identifier
    
# Subset based on gene IDs
sub <- gff %>% filter(ID %in% ids) %>% select(seqname, source, feature, start, end, score, strand, frame, attributes)
not <- gff %>% filter(ID %ni% ids)  %>% select(seqname, source, feature, start, end, score, strand, frame, attributes)

# Create Description
NCBIamr$description <- paste("Additional_database=NDARO;NDARO:Gene_Name=", NCBIamr$Gene.symbol, ";",
                     "NDARO:Gene_Product=", NCBIamr$Sequence.name, ";", "NDARO:Resistance_Category=",
                     NCBIamr$Element.type, ";", "NDARO:Resistance_Target=", NCBIamr$Class, ";",
                     "NDARO:Method=", NCBIamr$Method, ";", "NDARO:Closest_Sequence=", NCBIamr$Name.of.closest.sequence, sep = "")
NCBIamr$description <- gsub(" ", "_", NCBIamr$description)
    
## Add New Source
s <- sub$source
sn <- opt$database
snew <- paste(s, sn, sep = ",")
sub$source <- snew

## Add New Feature
f <- sub$feature
fn <- opt$type
fnew <- paste(f, fn, sep = ",")
sub$feature <- fnew

## attributes
sub$ID <- getAttributeField(as.character(sub$attributes), "ID", ";")
sub <- merge.data.frame(sub, NCBIamr, by.x = "ID", 
                        by.y = "Protein.identifier", all = TRUE)
sub <- unite(sub, "attributes", c("attributes", "description"), sep = ";") %>%
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
                col.names = FALSE, row.names = FALSE)

} else {
  # Load GFF file
  gff <- gffRead(opt$gff)
  # Write output
  write.table(gff, file = opt$out, quote = FALSE, sep = "\t", 
              col.names = FALSE, row.names = FALSE)
}
