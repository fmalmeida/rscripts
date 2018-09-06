#!/usr/bin/Rscript

# Load Libraries
suppressMessages(library(DataCombine))
suppressMessages(library(ballgown))
suppressMessages(library(optparse))

# Set parameters
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="output", 
              help="output file prefix name [default= %default]", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="feature type [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

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

# Load gff file
gff <- gffRead(opt$input)
output_file <- opt$out

## Create file specific for CARD database - Since it is the one that has the better described resistance genes

if (is.null(opt$type)) {
  ### Create fields - Prokka
  gff$Prokka_ID <- getAttributeField(gff$attributes, "ID", ";")
  gff$Prokka_gene <- getAttributeField(gff$attributes, "gene", ";")
  gff$Prokka_geneFamily <- substr(gff$Prokka_gene, 1, 3)
  gff$Prokka_name <- getAttributeField(gff$attributes, "Name", ";")
  gff$Prokka_product <- getAttributeField(gff$attributes, "product", ";")
  gff$Prokka_inference <- getAttributeField(gff$attributes, "inference", ";")
  gff$Domain <- getAttributeField(gff$attributes, "protein_motif", ";")
  gff$Additional_DB <- getAttributeField(gff$attributes, "Additional_database", ";")
  
  col = c("seqname", "start", "end", "feature", "source", "Prokka_ID", "Prokka_gene", "Prokka_geneFamily", 
          "Prokka_name", "Prokka_product", "Prokka_inference", "Domain", "Additional_DB")
  
  table <- gff[, col]
  out <- paste0(output_file, ".tsv", sep = "")
  write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else if (opt$type == "resistance") {
  
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
gff$Prokka_gene <- getAttributeField(gff$attributes, "gene", ";")
gff$Prokka_geneFamily <- substr(gff$Prokka_gene, 1, 3)
gff$Prokka_name <- getAttributeField(gff$attributes, "Name", ";")
gff$Prokka_product <- getAttributeField(gff$attributes, "product", ";")
gff$Prokka_inference <- getAttributeField(gff$attributes, "inference", ";")
gff$Domain <- getAttributeField(gff$attributes, "protein_motif", ";")
gff$Additional_DB <- getAttributeField(gff$attributes, "Additional_database", ";")

col = c("seqname", "start", "end", "feature", "source", "Prokka_ID", "Prokka_gene", "Prokka_geneFamily", 
        "Prokka_name", "Prokka_product", "Prokka_inference", "Domain", "Additional_DB")

table <- gff[, col]
out <- paste0(output_file, "_", opt$type, ".tsv", sep = "")
write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)}