#!/usr/bin/Rscript
# Setting Help
'usage: addBlast2Gff.R [--input=<file> --gff=<file> --out=<chr> --database=<chr> --type=<chr> --scoverage=<int>]

options:
  -i, --input=<file>    Tabular Blast to be added to GFF
  -g, --gff=<file>      GFF file to add Blast hits into
  -o, --out=<chr>       Output file name [default: out.gff]
  -d, --database=<chr>  Name of databased which Blast came from
  -t, --type=<chr>      Type of feature blasted. Ex: resistance
  -c, --scoverage=<int> Minimum subject coverage to keep' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)

if (is.null(opt$input)){
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

# Load libraries
suppressMessages(library(ballgown))
suppressMessages(library(DataCombine))
suppressMessages(library(dplyr))

# Function used to remove redundancy
reduce_row = function(i) {
  d <- unlist(strsplit(i, split=","))
  paste(unique(d), collapse = ',') 
}
# Function to get Attribute Fields
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
# Operator to discard patterns found
'%ni%' <- Negate('%in%')

# Check if file is empty

if (file.info(opt$input)$size > 0 ) {
# Load blast tabular file
blastHeader <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
                 "qend", "sstart", "send", "slen", "evalue", "bitscore", "stitle")

blastFile <- read.delim(opt$input, header = FALSE)
colnames(blastFile) <- blastHeader

# Filter blast based on subject coverage
if (!is.null(opt$scoverage)) {
blastFile$scov <- (blastFile$length / blastFile$slen) * 100
blastFile <- dplyr::filter(blastFile, scov >= as.integer(opt$scoverage))
}

if (nrow(blastFile) > 0) {
  
# Remove duplicates based on bitscore
blastFile <- blastFile[order(blastFile$qseqid, -abs(blastFile$bitscore) ), ]
blastFile <- dplyr::filter(blastFile, bitscore >= 500)
blastFile <-blastFile[ !duplicated(blastFile$qseqid), ]
blastFile <- blastFile[order(blastFile$qseqid),]

# Create GFF Attribute Entry
att <- paste("Additional_database=", opt$database, ";", opt$database, "_ID=", 
             blastFile$sseqid, ";", opt$database, "_Target=", blastFile$stitle, sep = "")

# Get gene names
ids <- blastFile$qseqid

# Load GFF file
gff <- gffRead(opt$gff)

# Create a column in gff with ids
gff$ID <- getAttributeField(gff$attributes, "ID", ";")

# Subset based on gene IDs
sub <- gff %>% filter(ID %in% ids) %>% select(seqname, source, feature, start, end, score, strand, frame, attributes)
not <- gff %>% filter(ID %ni% ids) %>% select(seqname, source, feature, start, end, score, strand, frame, attributes)

# Change fields values
## source
s <- sub$source
sn <- opt$database
snew <- paste(s, sn, sep = ",")
sub$source <- snew

## feature
f <- sub$feature
fn <- opt$type
fnew <- paste(f, fn, sep = ",")
sub$feature <- fnew

## attributes
a <- sub$attributes
an <- att
anew <- paste(a, an, sep = ";")
sub$attributes <- anew

# Merge files
merged_df <- merge.data.frame(sub, not, all = TRUE)
feat <- merged_df$feature
merged_df$feature <- sapply(feat, reduce_row)
source <- merged_df$source
merged_df$source <- sapply(source, reduce_row)
merged_df <- merged_df[order(merged_df$seqname, merged_df$start),]

# Write output
write.table(merged_df, file = opt$out, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE)
} else {
  # Load GFF file
  gff <- gffRead(opt$gff)
  # Write output
  write.table(gff, file = opt$out, quote = FALSE, sep = "\t", 
              col.names = FALSE, row.names = FALSE)
}} else {
  # Load GFF file
  gff <- gffRead(opt$gff)
  # Write output
  write.table(gff, file = opt$out, quote = FALSE, sep = "\t", 
              col.names = FALSE, row.names = FALSE)
}