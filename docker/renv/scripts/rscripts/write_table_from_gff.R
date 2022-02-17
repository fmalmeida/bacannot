#!/usr/bin/Rscript
# Setting help
'usage: write_table_from_gff.R [--input=<file> --out=<chr> --type=<chr>]

options:
  -i, --input=<file>    GFF file name
  -o, --out=<chr>       Output prefix file name [default: out]
  -t, --type=<chr>      Feature type to subset and write table from [default: resistance]' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)
if (is.null(opt$input)){
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Load Libraries
suppressMessages(library(DataCombine))
suppressMessages(library(ballgown))
suppressMessages(library(stringr))
  
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

getAdditionalProducts <- function (vector) {
  s=strsplit(vector, split = ";", fixed = TRUE)
  sapply(s, function(x) {
    v = str_subset(x, pattern="_target")
    w = strsplit(v, split = "=", fixed = TRUE)
    d = sapply(sapply(w, "[", 1), 
               function(x) {
                 strsplit(x, split = "_", fixed = TRUE)
                 })
    d = sapply(d, "[", 1)
    j = sapply(w, "[", 2)
    i = paste(d, j, sep = ":")
    if (length(i) > 0) {
      rv = paste(unique(i), collapse = ",") } else {rv = as.character(NA)}
  })
}

# Check if file is empty
if (file.info(opt$input)$size > 0) {
# Load gff file
gff <- gffRead(opt$input)
gff$attributes <- gsub(x = gff$attributes, pattern = ",id", replacement = ";id")
output_file <- opt$out

if (length(opt$type) && opt$type != "card") {
  ### Create fields - Prokka
  gff$Prokka_ID <- getAttributeField(gff$attributes, "id", ";")
  gff$Prokka_product <- getAttributeField(gff$attributes, "product", ";")
  gff$Prokka_inference <- getAttributeField(gff$attributes, "inference", ";")
  gff$Domain <- getAttributeField(gff$attributes, "protein_motif", ";")
  gff$Additional_DB <- getAttributeField(gff$attributes, "additional_database", ";")
  gff$Additional_product <- getAdditionalProducts(gff$attributes)
  
  ### Give columns a name
  col = c("seqname", "Prokka_ID", "start", "end", "feature", "source", "Additional_DB",
          "Prokka_product", "Additional_product", "Prokka_inference", "Domain")
  
  ### Write document
  table <- gff[, col]
  out <- paste0(output_file, "_", opt$type, ".tsv", sep = "")
  write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else if (opt$type == "card") {

  ### Create CARD specific summary table, 
  ### since it has great information about antimicrobial genes.
  gff$ARO_Accession <- getAttributeField(gff$attributes, "aro", ";")
  gff$Gene_Family <- getAttributeField(gff$attributes, "gene_family", ";")
  gff$Name <- getAttributeField(gff$attributes, "db_name", ";")
  gff$Drug_Class <- getAttributeField(gff$attributes, "drug_class", ";")
  gff$Resistance_Mechanism <- getAttributeField(gff$attributes, "resistance_mechanism", ";")
  gff$CVTerm <- getAttributeField(gff$attributes, "cvterm", ";")
  gff$Domain <- getAttributeField(gff$attributes, "protein_motif", ";")
  gff$Prokka_product <- getAttributeField(gff$attributes, "product", ";")
  gff$Prokka_ID <- getAttributeField(gff$attributes, "id", ";")
  gff$Prokka_inference <- getAttributeField(gff$attributes, "inference", ";")
  gff$Additional_DB <- getAttributeField(gff$attributes, "additional_database", ";")
  gff$Additional_product <- getAdditionalProducts(gff$attributes)

  #### Give columns a name
  col = c("seqname", "start", "end", "feature", "source", "ARO_Accession", 
        "Gene_Family", "Name", "Drug_Class", "Resistance_Mechanism", "CVTerm",
        "Domain", "Prokka_product", "Prokka_ID", "Prokka_inference", "Additional_DB", "Additional_product")
  
  ### Write document
  table <- gff[, col]
  out <- paste0(output_file, "_", opt$type, ".tsv", sep = "")
  write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
  ### Create non-specific file
  ### Create fields - Prokka
  gff$Prokka_ID <- getAttributeField(gff$attributes, "id", ";")
  gff$Prokka_product <- getAttributeField(gff$attributes, "product", ";")
  gff$Prokka_inference <- getAttributeField(gff$attributes, "inference", ";")
  gff$Annotated_domain <- getAttributeField(gff$attributes, "protein_motif", ";")

  ### Give columns a name
  col = c("seqname", "start", "end", "feature", "source", "Prokka_ID", 
        "Prokka_product", "Prokka_inference", "Annotated_domain")


  ### Write document
  table <- gff[, col]
  out <- paste0(output_file, "_non_specific.tsv", sep = "")
  write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}} else {
  opt <- options(show.error.messages=FALSE)
  on.exit(options(opt))
  stop()
}
