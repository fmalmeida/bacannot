#!/usr/bin/Rscript
# Setting Help
'usage: addBlast2Gff.R [--input=<file> --gff=<file> --out=<chr> --database=<chr> --type=<chr>]

options:
  -i, --input=<file>    Tabular Blast to be added to GFF
  -g, --gff=<file>      GFF file to add Blast hits into
  -o, --out=<chr>       Output file name [default: out.gff]
  -d, --database=<chr>  Name of databased which Blast came from
  -t, --type=<chr>      Type of feature blasted. Ex: resistance' -> doc

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
suppressMessages(library(tidyr))
suppressMessages(library(stringr))

#################
### Functions ###
#################

# Function used to remove redundancy
reduce_row = function(i) {
  d <- unlist(strsplit(i, split=","))
  paste(unique(d), collapse = ',')
}

# Parse Blast Titles
subj_title = function(x, db) {
  desc = strsplit(x, "~~~", fixed=TRUE)
  
  if (db == "PHAST" | db == "Victors") {
    text <- paste("Additional_database=", desc[[1]][1], ";", desc[[1]][1], "_Target=",
                  desc[[1]][4], ";", desc[[1]][1], "_Product=", desc[[1]][2], sep = "")
  } else if (db == "VFDB") {
    text <- paste("Additional_database=", desc[[1]][1], ";", desc[[1]][1], "_Target=",
                  desc[[1]][3], ";", desc[[1]][1], "_Product=", desc[[1]][4], sep = "")
  } else {
    text <- paste("Additional_database=", desc[[1]][1], ";", desc[[1]][1], "_Target=",
                  desc[[1]][2], ";", desc[[1]][1], "_Product=", desc[[1]][4], sep = "")
  }
  
  return(text)
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

#############
### BEGIN ###
#############

# Check if file is empty

if (file.info(opt$input)$size > 0 ) {
  
  blastFile <- read.delim(opt$input, header = TRUE)

  if (nrow(blastFile) > 0) {
    
    # Sort entries
    blastFile <- blastFile[order(blastFile$qseqid),]
    
    # Remove whitespaces for GFF
    blastFile$sseqid <- gsub(" ", "_", x = blastFile$sseqid)

    # Create GFF Attribute Entry
    blastFile$NEW_attributes <- sapply(blastFile$sseqid, subj_title, db=opt$database)

    # Get gene names
    ids <- blastFile$qseqid

    # Load GFF file for merge
    gff <- gffRead(opt$gff)

    # Create a column in gff with ids
    gff$ID <- getAttributeField(gff$attributes, "ID", ";")

    # Subset based on gene IDs
    ## Lines with our IDs
    sub <- gff %>% 
      filter(ID %in% ids) %>% 
      select(seqname, source, feature, start, end, score, strand, frame, attributes, ID)
    ## Lines without our IDs
    not <- gff %>% 
      filter(ID %ni% ids) %>% 
      select(seqname, source, feature, start, end, score, strand, frame, attributes)

    # Change fields values
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
    sub <- merge.data.frame(sub, blastFile, by.x = "ID",
                            by.y = "qseqid", all = TRUE)
    sub <- unite(sub, "attributes", c("attributes", "NEW_attributes"), sep = ";") %>%
      select(seqname, source, feature, start, end, score, strand, frame, attributes)

    # Merge files
    merged_df <- merge.data.frame(sub, not, all = TRUE)
    feat <- merged_df$feature
    merged_df$feature <- sapply(feat, reduce_row)
    source <- merged_df$source
    merged_df$source <- sapply(source, reduce_row)
    merged_df <- merged_df[str_order(merged_df$attributes, numeric = TRUE), ]

    # Write output
    write.table(merged_df, file = opt$out, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
    
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
