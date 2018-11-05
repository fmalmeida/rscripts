#!/usr/bin/Rscript
# Setting help
'usage: write_table_from_gff.R [--input=<file> --out=<chr> --type=<chr>]

options:
  -i, --input=<file>    GFF file name
  -o, --out=<chr>       Output prefix file name [default: out]
  -t, --type=<chr>      Feature type to subset and write table from [default: resistance]' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)
if (is.null(opt$input)){
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Load Libraries
suppressMessages(library(DataCombine))
suppressMessages(library(ballgown))

# Function used to get values from attributes colum
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

# Check if file is empty
if (file.info(opt$input)$size > 0) {
# Load gff file
gff <- gffRead(opt$input)
output_file <- opt$out

## Create non-specific file
  ### Create fields - Prokka
  gff$Prokka_ID <- getAttributeField(gff$attributes, "ID", ";")
  gff$Prokka_product <- getAttributeField(gff$attributes, "product", ";")
  gff$Prokka_inference <- getAttributeField(gff$attributes, "inference", ";")
  gff$Domain <- getAttributeField(gff$attributes, "protein_motif", ";")
  gff$Additional_DB <- getAttributeField(gff$attributes, "Additional_database", ";")
  gff$Additional_product <- getAttributeField(gff$attributes, "", ";")
  
  col = c("seqname", "start", "end", "feature", "source", "Prokka_ID", 
          "Prokka_product", "Prokka_inference", "Domain", "Additional_DB")
  
  table <- gff[, col]
  out <- paste0(output_file, ".tsv", sep = "")
  write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else if (opt$type == "CARD") {

# Create CARD specific summary table, since it has great information about antimicrobial genes.
    
  gff$ARO_Accession <- getAttributeField(gff$attributes, "ARO", ";")
  gff$Gene_Family <- getAttributeField(gff$attributes, "Gene_Family", ";")
  gff$Name <- getAttributeField(gff$attributes, "DB_Name", ";")
  gff$Drug_Class <- getAttributeField(gff$attributes, "Drug_Class", ";")
  gff$Resistance_Mechanism <- getAttributeField(gff$attributes, "Resistance_Mechanism", ";")
  gff$Domain <- getAttributeField(gff$attributes, "protein_motif", ";")
  gff$Prokka_product <- getAttributeField(gff$attributes, "product", ";")
  gff$Prokka_inference <- getAttributeField(gff$attributes, "inference", ";")

#### Write table
col = c("seqname", "start", "end", "feature", "source", "ARO_Accession", 
        "Gene_Family", "Name", "Drug_Class", "Resistance_Mechanism", 
        "Domain", "Prokka_product", "Prokka_inference")
table <- gff[, col]
out <- paste0(output_file, "_", opt$type, ".tsv", sep = "")
write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
### Create fields - Prokka
gff$Prokka_ID <- getAttributeField(gff$attributes, "ID", ";")
gff$Prokka_product <- getAttributeField(gff$attributes, "product", ";")
gff$Prokka_inference <- getAttributeField(gff$attributes, "inference", ";")
gff$Domain <- getAttributeField(gff$attributes, "protein_motif", ";")
gff$Additional_DB <- getAttributeField(gff$attributes, "Additional_database", ";")

col = c("seqname", "start", "end", "feature", "source", "Prokka_ID", 
        "Prokka_product", "Prokka_inference", "Domain", "Additional_DB")

table <- gff[, col]
out <- paste0(output_file, "_", opt$type, ".tsv", sep = "")
write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
  opt <- options(show.error.messages=FALSE)
  on.exit(options(opt))
  stop()
}