#!/usr/bin/Rscript
# Setting Help
'usage: addNCBIamr2Gff.R [--gff=<file> --out=<chr> --database=<chr> --type=<chr>]

options:
-g, --gff=<file>      GFF file to add NCBI AMR Annotations into
-o, --out=<chr>       Output file name [default: out.gff]
-t, --type=<chr>      Type of feature. Ex: resistance
-d, --database=<chr>  Name of databased which Blast came from' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)

if (is.null(opt$gff)){
  stop("At least one argument must be supplied (gff file)\n", call.=FALSE)
}

# Load libraries
suppressMessages(library(ballgown))
suppressMessages(library(DataCombine))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

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

# Load GFF File
gff <- gffRead(opt$gff)
gff$ID <- getAttributeField(gff$attributes, "id", ";")

# Filter Entries that are annotated from NCBI AMR hmm
NCBIamr <- filter(gff, str_detect(attributes, "ncbifam-amr"))

if (is.null(NCBIamr) == FALSE & dim(NCBIamr)[1] != 0) {
  
# Get its ids
NCBIamr$ID <- getAttributeField(NCBIamr$attributes, "id", ";")
ids <- NCBIamr$ID
    
# Subset based on gene IDs
sub <- gff %>% filter(ID %in% ids) %>% select(seqname, source, feature, start, end, score, strand, frame, attributes)
not <- gff %>% filter(ID %ni% ids)  %>% select(seqname, source, feature, start, end, score, strand, frame, attributes)
    
## Add New Source
s <- sub$source
sn <- opt$database
snew <- paste(s, sn, sep = ",")
sub$source <- snew

## Add New Feature
f <- sub$feature
fn <- opt$type
fnew <- paste(f, fn, sep = ",")
sub$feature <- fnew

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
}
