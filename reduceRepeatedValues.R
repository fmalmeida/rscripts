#!/usr/bin/Rscript
# Setting help
'usage: reduceRepeatedValues.R [--input=<file> --out=<chr>]

options:
  -i, --input=<file>    GFF file
  -o, --out=<chr>       Output file name [default: out.gff]' -> doc

# Load libraries
suppressMessages(library(ballgown))
suppressMessages(library(DataCombine))
suppressMessages(library(docopt))

# Parse parameters
opt <- docopt(doc)

# Function used to remove redundancy
reduce_row = function(i) {
  d <- unlist(strsplit(i, split=","))
  paste(unique(d), collapse = ',') 
}

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

# Load GFF file
gff <- gffRead(opt$input)

## Remove repeated values
feature <- gff$feature
gff$feature <- sapply(feature, reduce_row)
source <- gff$source
gff$source <- sapply(source, reduce_row)
reduced_df <- gff[order(gff$seqname, gff$start),]

# Write output
write.table(reduced_df, file = opt$out, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE)