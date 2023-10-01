#!/usr/bin/Rscript
doc <- 'usage: gff2sql.R [--input=<file> --out=<chr> --fasta=<file> --nucleotide=<file> --aminoacid=<file>]

options:
  -i, --input=<file>         GFF file to transform in SQL
  -o, --out=<chr>            SQL database name to output [default: out.sql]
  -n, --nucleotide=<file>    Takes in the nucleotide FASTA.
  -a, --aminoacid=<file>     Takes in the protein FASTA
  -f, --fasta=<file>         Takes in the genome FASTA'

# Loading required packages
suppressMessages(library("docopt"))
suppressMessages(library(RSQLite))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(DataCombine))
suppressMessages(library(Biostrings))

# Parse help
opt <- docopt(doc)

# Useful functions
## Query the 9th column
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

## Add table to SQL db
addTable <- function (con, sql, input) {
  ## Open db
  suppressWarnings(dbBegin(con))

  ## Send rule
  res <- suppressWarnings(dbSendQuery(con, sql))

  ## Insert data based on rule
  suppressWarnings(dbBind(res, input))
  suppressWarnings(dbFetch(res))
  suppressWarnings(dbClearResult(res))

  ## Close db
  suppressWarnings(dbCommit(con))
}

# Loading SQL database driver
drv <- dbDriver("SQLite")
print(opt$out)
con <- dbConnect(drv, dbname=opt$out)

#####################################
### First STEP load GENOME to sql ###
#####################################
fastaFile <- readDNAStringSet(opt$fasta)
seq_name = names(fastaFile)
#sequence = paste(fastaFile)
sequence_len = sapply(fastaFile, function(x) {
  length(x)[[1]]
})
genome <- data.frame(seq_name, sequence_len)
names(genome) <- c("Contig", "Length")

# Create SQL table for the genome sequence
suppressWarnings(dbGetQuery(con, "CREATE Table Genome (Contig TEXT, Length TEXT)"))
# Create sql rule
sql <- "INSERT INTO Genome VALUES ($Contig, $Length)"
# Add to SQL db
addTable(con, sql, genome)

###################################
### Second STEP load GFF to sql ###
###################################

# Loading GFF file
gff <- read.delim(opt$input, header = FALSE, stringsAsFactors = FALSE)
# Give data a header
names(gff) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
# Get IDs
gff$ID <- getAttributeField(as.character(gff$attributes), "ID", ";")
# Reorder columns
gff <- gff %>% select(chr, source, ID, feature, start, end, score, strand, frame, attributes)
# Create SQL table to store GFF data
suppressWarnings(dbGetQuery(con, "CREATE Table GFF (Contig TEXT, Source TEXT, ID TEXT, Feature TEXT,
           Start INTEGER, End INTEGER, Score INTEGER, Strand TEXT,
           Frame INTEGER, Attributes TEXT)"))
# Create sql rule
sql <- "INSERT INTO GFF VALUES ($chr, $source, $ID, $feature,
$start, $end, $score, $strand, $frame, $attributes)"
# Add to SQL db
addTable(con, sql, gff)

##############################################
### Third STEP load gene nucl fasta to sql ###
##############################################

## Loading Protein fasta
genes <- readAAStringSet(opt$aminoacid)
gene_ids <- sapply(names(genes), function(x) {
  unlist(strsplit(as.character(x), " "))[1]
})
gene_desc <- sapply(names(genes), function(x) {
  paste0(unlist(strsplit(as.character(x), " "))[-1], collapse = " ")
})
sequences = paste(genes)
genes_aa <- data.frame(gene_ids, gene_desc, sequences)
names(genes_aa) <- c("ID", "Description", "Sequence")
## Create SQL table to store Protein FASTA
suppressWarnings(dbGetQuery(con, "CREATE Table ProteinFasta (ID TEXT, Description TEXT, Sequence TEXT)"))
## Create sql rule
sql <- "INSERT INTO ProteinFasta VALUES ($ID, $Description, $Sequence)"
# Add to SQL db
addTable(con, sql, genes_aa)

## Loading Nucleotide fasta
genes <- readDNAStringSet(opt$nucleotide)
gene_ids <- sapply(names(genes), function(x) {
  unlist(strsplit(as.character(x), " "))[1]
})
gene_desc <- sapply(names(genes), function(x) {
  paste0(unlist(strsplit(as.character(x), " "))[-1], collapse = " ")
})
sequences = paste(genes)
genes_ncl <- data.frame(gene_ids, gene_desc, sequences)
names(genes_ncl) <- c("ID", "Description", "Sequence")
## Create SQL table to store Protein FASTA
suppressWarnings(dbGetQuery(con, "CREATE Table NucleotideFasta (ID TEXT, Description TEXT, Sequence TEXT)"))
## Create sql rule
sql <- "INSERT INTO NucleotideFasta VALUES ($ID, $Description, $Sequence)"
# Add to SQL db
addTable(con, sql, genes_ncl)
