#!/usr/bin/Rscript

# Load CARD entries index
cat_index <- read.table("/work/indexes/aro_categories_index.csv", header = TRUE, sep = "\t")
cat <- read.table("/work/indexes/aro_categories.csv", header = TRUE, sep = "\t")
index <- read.table("/work/indexes/aro_index.csv", header = TRUE, sep = "\t", fill = TRUE)

# Load library
suppressMessages(library(DataCombine))
suppressMessages(library(ballgown))

# Setting parameters
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="blast tab output", metavar="character"),
  make_option(c("-g", "--gff"), type="character", default = NULL,
              help = "gff file to merge", metavar = "character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-d", "--database"), type="character", default=NULL,
              help="database annotated from", metavar="character"),
  make_option(c("-t", "--type", type="character"), default = NULL,
              help = "feature type", metavar = "character"),
  make_option(c("-p", "--pident"), default=90, metavar = "integer",
              help = "% identitty to filter blast [default= %default]", type = "integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

# Merge indexes to create a full index with all values
merged <- merge.data.frame(cat_index, index, by.x = "Protein.Accession",
                           by.y = "Protein.Accession", all = TRUE)
card_indexes <- merge.data.frame(merged, cat, by.x = "ARO.Accession", 
                                 by.y = "ARO.Accession", all = TRUE)

# Load Blast
blastFile <- read.table(opt$input, sep = "\t")
blastHeader <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(blastFile) <- blastHeader

# Remove duplicates based on bitscore
blastFile <- blastFile[order(blastFile$qseqid, -abs(blastFile$bitscore) ), ]
blastFile <-blastFile[ !duplicated(blastFile$qseqid), ]
blastFile <- blastFile[order(blastFile$qseqid),]

## Filter per %Identity
blast_filtered <- subset(blastFile, pident >= opt$pident)
ssids <- as.vector(blast_filtered$sseqid)
aroID <- sapply(strsplit(ssids, "\\|"), `[`, 3)
blast_filtered$ARO <- aroID

#Load gff
gff <- gffRead(opt$gff)

##Subset Card indexes
card_subset <- grepl.sub(card_indexes, pattern = aroID, Var = "ARO.Accession")

#Get desired fields
description <- paste("Additional_Database=", opt$database, ";", 
                     "ARO=", card_subset$ARO.Accession, ";", "Gene_Family=", 
                     card_subset$AMR.Gene.Family, ";", "Drug_Class=", card_subset$Drug.Class, 
                     ";", "Resistance_Mechanism=", card_subset$Resistance.Mechanism, ";", 
                     "DB_Name=", card_subset$Model.Name, sep = "")

card_subset$attributes <- description

# Concatenate new attributes values
blast_filtered <- merge.data.frame(blast_filtered, card_subset, by.x = "ARO", 
                                   by.y = "ARO.Accession", all = TRUE)

blast_filtered <- blast_filtered[order(blastFile$qseqid),]

# Get gene names from blast hits
prokka_ids <- blast_filtered$qseqid

# Subset GFF - based on genes that have a hit
sub <- grepl.sub(gff, pattern = prokka_ids, Var = "attributes")
not <- grepl.sub(gff, pattern = prokka_ids, Var = "attributes", keep.found = FALSE)

#Change fields - Add database source and feature type
##source
s <- sub$source
sn <- opt$database
snew <- paste(s, sn, sep = ",")
sub$source <- snew

##feature
f <- sub$feature
fn <- opt$type
fnew <- paste(f, fn, sep = ",")
sub$feature <- fnew

##attributes
blast_filtered <- blast_filtered[order(blast_filtered$qseqid),]
a <- sub$attributes
an <- blast_filtered$attributes
anew <- paste(a, an, sep = ";")
sub$attributes <- anew

##Merging final GFF
reduce_row = function(i) {
  d <- unlist(strsplit(i, split=","))
  paste(unique(d), collapse = ',') 
}

#Merge files
merged_df <- merge.data.frame(sub, not, all = TRUE)
feat <- merged_df$feature
merged_df$feature <- sapply(feat, reduce_row)
source <- merged_df$source
merged_df$source <- sapply(source, reduce_row)
merged_df <- merged_df[order(merged_df$seqname, merged_df$start),]

# Write output
write.table(merged_df, file = opt$out, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE, append = FALSE)

#Clear workspace
rm(list=ls())