#!/usr/bin/Rscript

suppressMessages(library(ballgown))
suppressMessages(library(DataCombine))

# Set function
getmotif <- function (x, field, attrsep = ";") { 
  s = strsplit(x, split = attrsep, fixed = TRUE) 
  sapply(s, function(atts) { 
    a = strsplit(atts, split = ":", fixed = TRUE) 
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

gff <- gffRead(opt$input)

# Victors
sub <- grepl.sub(gff, pattern = "*victors*", Var = "attributes")
not <- grepl.sub(gff, pattern = "*victors*", Var = "attributes", keep.found = FALSE)
# source
s <- sub$source
sn <- "victors"
snew <- paste(s, sn, sep = ",")
sub$source <- snew
# features
f <- sub$feature
fn <- "virulence"
fnew <- paste(f, fn, sep = ",")
sub$feature <- fnew
# merge
merged_df <- merge.data.frame(sub, not, all = TRUE)
feat <- merged_df$feature
merged_df$feature <- sapply(feat, reduce_row)
source <- merged_df$source
merged_df$source <- sapply(source, reduce_row)
merged_df <- merged_df[order(merged_df$seqname, merged_df$start),]

# VFDB
sub <- grepl.sub(merged_df, pattern = "*VFDB*", Var = "attributes")
not <- grepl.sub(merged_df, pattern = "*VFDB*", Var = "attributes", keep.found = FALSE)
# source
s <- sub$source
sn <- "VFDB"
snew <- paste(s, sn, sep = ",")
sub$source <- snew
# features
f <- sub$feature
fn <- "virulence"
fnew <- paste(f, fn, sep = ",")
sub$feature <- fnew
# merge
merged_df <- merge.data.frame(sub, not, all = TRUE)
feat <- merged_df$feature
merged_df$feature <- sapply(feat, reduce_row)
source <- merged_df$source
merged_df$source <- sapply(source, reduce_row)
merged_df <- merged_df[order(merged_df$seqname, merged_df$start),]

# ICEberg
sub <- grepl.sub(merged_df, pattern = "*ICEberg*", Var = "attributes")
not <- grepl.sub(merged_df, pattern = "*ICEberg*", Var = "attributes", keep.found = FALSE)
# source
s <- sub$source
sn <- "ICEberg"
snew <- paste(s, sn, sep = ",")
sub$source <- snew
# features
f <- sub$feature
fn <- "ICE"
fnew <- paste(f, fn, sep = ",")
sub$feature <- fnew
# merge
merged_df <- merge.data.frame(sub, not, all = TRUE)
feat <- merged_df$feature
merged_df$feature <- sapply(feat, reduce_row)
source <- merged_df$source
merged_df$source <- sapply(source, reduce_row)
merged_df <- merged_df[order(merged_df$seqname, merged_df$start),]

# PHAGE
sub <- grepl.sub(merged_df, pattern = "*prophage*", Var = "attributes")
not <- grepl.sub(merged_df, pattern = "*prophage*", Var = "attributes", keep.found = FALSE)
# source
s <- sub$source
sn <- "phast"
snew <- paste(s, sn, sep = ",")
sub$source <- snew
# features
f <- sub$feature
fn <- "PHAGE"
fnew <- paste(f, fn, sep = ",")
sub$feature <- fnew
# merge
merged_df <- merge.data.frame(sub, not, all = TRUE)
feat <- merged_df$feature
merged_df$feature <- sapply(feat, reduce_row)
source <- merged_df$source
merged_df$source <- sapply(source, reduce_row)
merged_df <- merged_df[order(merged_df$seqname, merged_df$start),]

# Resistance
sub <- grepl.sub(merged_df, pattern = "*resistance*", Var = "attributes")
not <- grepl.sub(merged_df, pattern = "*resistance*", Var = "attributes", keep.found = FALSE)
# features
f <- sub$feature
fn <- "resistance"
fnew <- paste(f, fn, sep = ",")
sub$feature <- fnew
# merge
merged_df <- merge.data.frame(sub, not, all = TRUE)
feat <- merged_df$feature
merged_df$feature <- sapply(feat, reduce_row)
source <- merged_df$source
merged_df$source <- sapply(source, reduce_row)
merged_df <- merged_df[order(merged_df$seqname, merged_df$start),]

# resfams
sub <- grepl.sub(merged_df, pattern = "*resfams*", Var = "attributes")
not <- grepl.sub(merged_df, pattern = "*resfams*", Var = "attributes", keep.found = FALSE)
# source
s <- sub$source
sn <- "Resfams"
snew <- paste(s, sn, sep = ",")
sub$source <- snew
# features
f <- sub$feature
fn <- "resistance"
fnew <- paste(f, fn, sep = ",")
sub$feature <- fnew
# merge
merged_df <- merge.data.frame(sub, not, all = TRUE)
feat <- merged_df$feature
merged_df$feature <- sapply(feat, reduce_row)
source <- merged_df$source
merged_df$source <- sapply(source, reduce_row)
merged_df <- merged_df[order(merged_df$seqname, merged_df$start),]

# Write output
write.table(merged_df, file = opt$out, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE)
