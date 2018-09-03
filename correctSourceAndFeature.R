#!/usr/bin/Rscript

suppressMessages(library(ballgown))
suppressMessages(library(DataCombine))

# Set function
# This function is used to extract the values of the fields stored in
# Attributes column of gff file.
getmotif <- function (x, field, attrsep = ";") { 
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

# Read gff file
gff <- gffRead(opt$input)
gff$attributes <- gsub(" ", "_", gff$attributes)

# First step is to subset those entries that have something related to one 
# of the pfam subsets

sub <- grepl.sub(gff, pattern = "_subset", Var = "attributes")
not <- grepl.sub(gff, pattern = "_subset", Var = "attributes", keep.found = FALSE)
sub$attributes <- gsub(",protein_motif:", ";protein_motif=", sub$attributes)
sub$attributes <- gsub(",", ";", sub$attributes)

# Secondly we need to store the previous value of the source column in 
# order to add to it the name of the pfam subset database

previous_source <- sub$source

# Then, we need to extract the name of the database.
motifs <- getmotif(sub$attributes, "protein_motif", ";")
motifs <- sapply(strsplit(motifs, split = "_", fixed = TRUE), "[", 1)

# Finally, we add these names to the previous sources
new_source <- paste(previous_source, motifs, sep = ",")
sub$source <- new_source

# Next, we need to write the respective feature types.
features <- gsub("victors|VFDB", "virulence", motifs)
features <- gsub("ICEberg", "ICE", features)

# Then, we add it to previously wrote features
previous_features <- sub$feature
new_features <- paste(previous_features, features, sep = ",")
sub$feature <- new_features

# Merge gff subsets
merged_df <- merge.data.frame(sub, not, all = TRUE)
feat <- merged_df$feature
merged_df$feature <- sapply(feat, reduce_row)
source <- merged_df$source
merged_df$source <- sapply(source, reduce_row)
merged_df <- merged_df[order(merged_df$seqname, merged_df$start),]

# Write output
write.table(merged_df, file = opt$out, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE)
