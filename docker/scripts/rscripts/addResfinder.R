#!/usr/bin/Rscript
# Setting Help
'usage: addResfinder.R [--txt=<file> --gff=<file> --out=<chr>]

options:
-g, --gff=<file>      GFF file to merge annotation
-t, --txt=<file>      Resfinder intersect file
-o, --out=<chr>       Output file name [default: out.gff]' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)

if (is.null(opt$gff)){
  stop("At least one argument must be supplied (gff file)\n", call.=FALSE)
}

if (is.null(opt$txt)){
  stop("At least one argument must be supplied (intersection file)\n", call.=FALSE)
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
  s = strsplit(as.character(x), split = attrsep, fixed = TRUE)
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

if (file.info(opt$txt)$size > 0) {
  
  # Load GFF file
  gff <- gffRead(opt$gff)
  
  # Create a column in the intersection file with ids
  gff$ID <- getAttributeField(gff$attributes, "ID", ";")
  
  # Load intersection file
  resfinder <- read.csv(opt$txt, header = F, sep = "\t")
  colnames(resfinder) <- c("seqname1", "source1", "feature1", "start1", "end1", "score1", "strand1", "frame1", "attributes1",
                     "seqname2", "source2", "feature2", "start2", "end2", "score2", "strand2", "frame2", "attributes2",
                     "len")
  
  # Create a column in the intersection file with ids
  resfinder$ID <- getAttributeField(resfinder$attributes2, "ID", ";")
  
  # save ids
  ids <- resfinder$ID

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
  sn <- "Resfinder"
  snew <- paste(s, sn, sep = ",")
  sub$source <- snew

  ## feature
  f <- sub$feature
  fn <- "Resistance"
  fnew <- paste(f, fn, sep = ",")
  sub$feature <- fnew

  ## attributes
  sub <- merge.data.frame(sub, resfinder, by = "ID", all = TRUE)
  sub$attributes1 <- gsub(pattern = "ID=", replacement = "Resfinder_ID=", x=sub$attributes1)
  sub <- unite(sub, "attributes", c("attributes", "attributes1"), sep = ";") %>%
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
  write.table(gff, file = opt$out, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
}
