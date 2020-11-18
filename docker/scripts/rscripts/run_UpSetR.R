#!/usr/bin/Rscript
# Setting help
'usage: run_UpSetR.R [--input=<file> --out=<chr> --mainColor=<chr> --setsColor=<chr>]

options:
  -i, --input=<file>    Input to Plot
  -o, --out=<chr>       Output name to plot  [default: out.png]
  --mainColor=<chr>     Sets Main Bars Color [default: gray23]
  --setsColor=<chr>     Sets horizontal bars Color [default: gray23]' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)
if (is.null(opt$input)){
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

# Load Library
suppressMessages(library(UpSetR))

# Get DataFrame
df <- read.delim(opt$input, header=TRUE)

# Run UpSetR
png(opt$out)
upset(df, order.by="freq", sets.bar.color=opt$setsColor, main.bar.color=opt$mainColor, keep.order = TRUE)
dev.off()
