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
suppressMessages(library(stringr))
  
# Function used to get values from attributes column. It was completely edited
# to keep all additional findings. I may have to implement on other scripts.
getAttributeField <- function (df, field, attrsep = ";") {
  s = strsplit(df, split = attrsep, fixed = TRUE)
  v <- sapply(s, function(x) {
    v = str_subset(x, pattern=field)
    y = strsplit(v, split = "=", fixed = TRUE)
    m = sapply(y, "[", 2)
    if (length(m) > 0) {
    rv = paste(unique(m), collapse = ",") } else { 
      rv = as.character(NA) 
    } 
  })
  return(v)
}

getAdditionalProducts <- function (vector) {
  s=strsplit(vector, split = ";", fixed = TRUE)
  sapply(s, function(x) {
    v = str_subset(x, pattern="_Target")
    w = strsplit(v, split = "=", fixed = TRUE)
    d = sapply(sapply(w, "[", 1), 
               function(x) {
                 strsplit(x, split = "_", fixed = TRUE)
                 })
    d = sapply(d, "[", 1)
    j = sapply(w, "[", 2)
    i = paste(d, j, sep = ":")
    if (length(i) > 0) {
      rv = paste(unique(i), collapse = ",") } else {rv = as.character(NA)}
  })
}

# Check if file is empty
if (file.info(opt$input)$size > 0) {
# Load gff file
gff <- gffRead(opt$input)
gff$attributes <- gsub(x = gff$attributes, pattern = ",ID", replacement = ";ID")
output_file <- opt$out

if (exists("opt$type") && opt$type != "CARD") {
  ### Create fields - Prokka
  gff$Prokka_ID <- getAttributeField(gff$attributes, "ID", ";")
  gff$Prokka_product <- getAttributeField(gff$attributes, "product", ";")
  gff$Prokka_inference <- getAttributeField(gff$attributes, "inference", ";")
  gff$Domain <- getAttributeField(gff$attributes, "protein_motif", ";")
  gff$Additional_DB <- getAttributeField(gff$attributes, "Additional_database", ";")
  gff$Additional_product <- getAdditionalProducts(gff$attributes)
  
  ### Give columns a name
  col = c("seqname", "Prokka_ID", "start", "end", "feature", "source", "Additional_DB",
          "Prokka_product", "Additional_product", "Prokka_inference", "Domain")
  
  ### Write document
  table <- gff[, col]
  out <- paste0(output_file, "_", opt$type, ".tsv", sep = "")
  write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else if (opt$type == "CARD") {

  ### Create CARD specific summary table, 
  ### since it has great information about antimicrobial genes.
  gff$ARO_Accession <- getAttributeField(gff$attributes, "ARO", ";")
  gff$Gene_Family <- getAttributeField(gff$attributes, "Gene_Family", ";")
  gff$Name <- getAttributeField(gff$attributes, "DB_Name", ";")
  gff$Drug_Class <- getAttributeField(gff$attributes, "Drug_Class", ";")
  gff$Resistance_Mechanism <- getAttributeField(gff$attributes, "Resistance_Mechanism", ";")
  gff$Domain <- getAttributeField(gff$attributes, "protein_motif", ";")
  gff$Prokka_product <- getAttributeField(gff$attributes, "product", ";")
  gff$Prokka_inference <- getAttributeField(gff$attributes, "inference", ";")
  gff$Additional_DB <- getAttributeField(gff$attributes, "Additional_database", ";")
  gff$Additional_product <- getAdditionalProducts(gff$attributes)

  #### Give columns a name
  col = c("seqname", "start", "end", "feature", "source", "ARO_Accession", 
        "Gene_Family", "Name", "Drug_Class", "Resistance_Mechanism", 
        "Domain", "Prokka_product", "Prokka_inference", "Additional_DB", "Additional_product")
  
  ### Write document
  table <- gff[, col]
  out <- paste0(output_file, "_", opt$type, ".tsv", sep = "")
  write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
  ### Create non-specific file
  ### Create fields - Prokka
  gff$Prokka_ID <- getAttributeField(gff$attributes, "ID", ";")
  gff$Prokka_product <- getAttributeField(gff$attributes, "product", ";")
  gff$Prokka_inference <- getAttributeField(gff$attributes, "inference", ";")
  gff$Annotated_domain <- getAttributeField(gff$attributes, "protein_motif", ";")

  ### Give columns a name
  col = c("seqname", "start", "end", "feature", "source", "Prokka_ID", 
        "Prokka_product", "Prokka_inference", "Annotated_domain")


  ### Write document
  table <- gff[, col]
  out <- paste0(output_file, ".tsv", sep = "")
  write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}} else {
  opt <- options(show.error.messages=FALSE)
  on.exit(options(opt))
  stop()
}