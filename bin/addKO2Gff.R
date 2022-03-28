#!/usr/bin/Rscript
# Setting Help
'usage: addKO2Gff.R [--input=<file> --gff=<file> --out=<chr> --database=<chr> --type=<chr> --scoverage=<int>]

options:
-i, --input=<file>    Tabular KOfamscan file to be added to GFF
-g, --gff=<file>      GFF file to add Blast hits into
-o, --out=<chr>       Output file name [default: out.gff]
-d, --database=<chr>  Name of databased which Blast came from' -> doc

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

# Check if file is empty

if (file.info(opt$input)$size > 0 ) {
  # Load KOfamscan tabular file
  KOfamHeader <- c("ID", "KO")
  
  KOfamFile <- read.delim(opt$input, header = FALSE)
  colnames(KOfamFile) <- KOfamHeader
  
  # Filter Only Proteins that have a KO
  KOfamFile <- dplyr::filter(KOfamFile, KO != "")
  
  if (nrow(KOfamFile) > 0) {
    
    # Create GFF Attribute Entry
    att <- paste("Additional_database=", opt$database, ";", "KO=", KOfamFile$KO, ";Method=KOfamscan", sep = "")
    
    # Get gene names
    ids <- KOfamFile$ID
    
    # Load GFF file
    gff <- gffRead(opt$gff)
    
    # Create a column in gff with ids
    gff$ID <- getAttributeField(gff$attributes, "ID", ";")
    
    # Subset based on gene IDs
    sub <- gff %>% filter(ID %in% ids) %>% select(seqname, source, feature, start, end, score, strand, frame, attributes)
    not <- gff %>% filter(ID %ni% ids)  %>% select(seqname, source, feature, start, end, score, strand, frame, attributes)
    
    ## attributes
    a <- sub$attributes
    an <- att
    anew <- paste(a, an, sep = ";")
    sub$attributes <- anew
    
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
  }} else {
    # Load GFF file
    gff <- gffRead(opt$gff)
    # Write output
    write.table(gff, file = opt$out, quote = FALSE, sep = "\t", 
                col.names = FALSE, row.names = FALSE)
  }