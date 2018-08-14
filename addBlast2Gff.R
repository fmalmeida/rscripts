#!/usr/bin/Rscript
suppressMessages(library(ballgown))
suppressMessages(library(DataCombine))
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
                help = "feature type", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

#Reduce row function
reduce_row = function(i) {
  d <- unlist(strsplit(i, split=","))
  paste(unique(d), collapse = ',') 
}

#Load blast result
blastHeader <- c("qseqid", "sseqid", "pident", "length", "mismatch",
"gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle")

blastFile <- read.table(opt$input, sep = "\t")
colnames(blastFile) <- blastHeader

att <- paste("Additional database=", opt$database, ";", "DB_ID=", blastFile$sseqid, ";DB_Target=",
             blastFile$stitle, sep = "")

ids <- blastFile$qseqid

#Load GFF file
gff <- gffRead(opt$gff)

#Subset
sub <- grepl.sub(gff, pattern = ids, Var = "attributes")
not <- grepl.sub(gff, pattern = ids, Var = "attributes", keep.found = FALSE)

#Change fields
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
a <- sub$attributes
an <- att
anew <- paste(a, an, sep = ";")
sub$attributes <- anew

#Merge files
merged_df <- merge.data.frame(sub, not, all = TRUE)
feat <- merged_df$feature
merged_df$feature <- sapply(feat, reduce_row)
source <- merged_df$source
merged_df$source <- sapply(source, reduce_row)
merged_df <- merged_df[order(merged_df$seqname, merged_df$start),]

# Write output
write.table(merged_df, file = opt$out, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE)

#Clear workspace
rm(list=ls())