#!/usr/bin/Rscript

# Load Library
suppressMessages(library(DataCombine))
suppressMessages(library(ballgown))

# optparse
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="output", 
              help="output file prefix name [default= %default]", metavar="character"),
  make_option(c("-p", "--pattern"), type = "character", default=NULL,
              help="subset pattern [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Create function
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

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
if (is.null(opt$pattern)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (pattern character).n", call.=FALSE)
}

# Load gff
gff_file <- opt$input
gff <- gffRead(gff_file)
output_file <- opt$out

# Subset (optional)
if ( !is.null(opt$pattern) ) {
  output_file <- opt$out 
  #sub_df <- grepl.sub(gff, pattern = opt$pattern, Var = "attributes")
  #card_df <- grepl.sub(gff, pattern = opt$pattern, Var = "feature")
  #merged_df <- merge.data.frame(sub_df, card_df, all = TRUE)
  merged_df <- grepl.sub(gff, pattern = opt$pattern, Var = "feature")
  merged_df$attributes <- gsub(",", ";", merged_df$attributes)
  sub_df <- merged_df
  #Create Fields
  sub_df$ID <- getAttributeField(sub_df$attributes, "ID", ";")
  sub_df$gene <- getAttributeField(sub_df$attributes, "gene", ";")
  sub_df$geneFamily <- substr(sub_df$gene, 1, 3)
  sub_df$name <- getAttributeField(sub_df$attributes, "Name", ";")
  sub_df$product <- getAttributeField(sub_df$attributes, "product", ";")
  # Write table
  col = c("seqname", "start", "end", "feature", "source", "ID", "gene", "geneFamily", "product")
  table <- sub_df[, col]
  out <- paste0(output_file, "_", opt$pattern, ".txt", sep = "")
  write.table(table, out, quote = FALSE, sep = "\t", 
              row.names = FALSE, col.names = TRUE)
}

# Output Name
output_file <- opt$out

# Create fields
gff$attributes <- gsub(",", ";", gff$attributes)
gff$ID <- getAttributeField(gff$attributes, "ID", ";")
gff$gene <- getAttributeField(gff$attributes, "gene", ";")
gff$geneFamily <- substr(gff$gene, 1, 3)
gff$name <- getAttributeField(gff$attributes, "Name", ";")
gff$product <- getAttributeField(gff$attributes, "product", ";")

# Write table
col = c("seqname", "start", "end", "feature", "source", "ID", "gene", "geneFamily", "product")
table <- gff[, col]
out <- paste0(output_file, ".txt", sep = "")
write.table(table, out, quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)
