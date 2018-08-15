#!/usr/bin/Rscript

suppressMessages(library(ballgown))
suppressMessages(library(DataCombine))

# Set function
reduce_row = function(i) {
  d <- unlist(strsplit(i, split=","))
  paste(unique(d), collapse = ',') 
}

# Setting parameters
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

#Load GFF
gff <- gffRead(opt$input)

## Reduce values
feature <- gff$feature
gff$feature <- sapply(feature, reduce_row)
source <- gff$source
gff$source <- sapply(source, reduce_row)
reduced_df <- gff[order(gff$seqname, gff$start),]

# Write output
write.table(reduced_df, file = opt$out, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE)
