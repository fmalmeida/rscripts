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
              help="output file prefix name [default= %default]", metavar="character")
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

# Load gff
gff_file <- opt$input
gff <- gffRead(gff_file)
output_file <- opt$out

# Genes annotated related to resistance
output_file <- opt$out
resistance_df <- grepl.sub(gff, pattern = "resistance", Var = "feature")

## Create file specific for CARD database - Since it is the one that has the better described genes

### CARD
card_df <- grepl.sub(resistance_df, pattern = "CARD", Var = "source")
card_df$ARO_Accession <- getAttributeField(card_df$attributes, "ARO", ";")
card_df$Gene_Family <- getAttributeField(card_df$attributes, "Gene_Family", ";")
card_df$Name <- getAttributeField(card_df$attributes, "DB_Name", ";")
card_df$Drug_Class <- getAttributeField(card_df$attributes, "Drug_Class", ";")
card_df$Resistance_Mechanism <- 
  getAttributeField(card_df$attributes, "Resistance_Mechanism", ";")

#### Write table
col = c("seqname", "start", "end", "feature", "source", "ARO_Accession", 
        "Gene_Family", "Name", "Drug_Class", "Resistance_Mechanism")
table <- card_df[, col]
out <- paste0(output_file, "_CARD_genes", ".tsv", sep = "")
write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

## Create a gene table that has all the databases or features together. With one specific column for 
## each database entry.

### Output Name
output_file <- opt$out

### Create fields - Prokka
gff$Prokka_ID <- getAttributeField(gff$attributes, "ID", ";")
gff$Prokka_gene <- getAttributeField(gff$attributes, "gene", ";")
gff$Prokka_geneFamily <- substr(gff$Prokka_gene, 1, 3)
gff$Prokka_name <- getAttributeField(gff$attributes, "Name", ";")
gff$Prokka_product <- getAttributeField(gff$attributes, "product", ";")
gff$Prokka_inference <- getAttributeField(gff$attributes, "inference", ";")

### Create fields - vfdb
gff$VFDB_ID <- getAttributeField(gff$attributes, "VFDB_ID", ";")
gff$VFDB_Target <- getAttributeField(gff$attributes, "VFDB_Target", ";")

### Create fields - victors
gff$Victors_ID <- getAttributeField(gff$attributes, "victors_ID", ";")
gff$Victors_Target <- getAttributeField(gff$attributes, "victors_Target", ";")

col = c("seqname", "start", "end", "feature", "source", "Prokka_ID", "Prokka_gene", "Prokka_geneFamily", 
        "Prokka_name", "Prokka_product", "Prokka_inference", "VFDB_ID", "VFDB_Target", "Victors_ID", "Victors_Target")
virulence <- grepl.sub(gff, "virulence", "feature")
table <- virulence[, col]
out <- paste0(output_file, "_virulence_gene_table", ".tsv", sep = "")
write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

### Create fields - CARD
gff$CARD_ARO_Accession <- getAttributeField(gff$attributes, "ARO", ";")
gff$CARD_Drug_Class <- getAttributeField(gff$attributes, "Drug_Class", ";")

### Create fields - Resfinder
gff$Resfinder_ID <- getAttributeField(gff$attributes, "resfinder_ID", ";")
gff$Resfinder_Target <- getAttributeField(gff$attributes, "resfinder_Target", ";")

col = c("seqname", "start", "end", "feature", "source", "Prokka_ID", "Prokka_gene", "Prokka_geneFamily", 
        "Prokka_name", "Prokka_product", "Prokka_inference", "CARD_ARO_Accession", "CARD_Drug_Class", 
        "Resfinder_ID", "Resfinder_Target")
resistance <- grepl.sub(gff, "resistance", "feature")
table <- resistance[, col]
out <- paste0(output_file, "_resistance_gene_table", ".tsv", sep = "")
write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

### Create fields - ICEberg
gff$ICEberg_ID <- getAttributeField(gff$attributes, "ICEberg_ID", ";")
gff$ICEberg_Target <- getAttributeField(gff$attributes, "ICEberg_Target", ";")

col = c("seqname", "start", "end", "feature", "source", "Prokka_ID", "Prokka_gene", "Prokka_geneFamily", 
        "Prokka_name", "Prokka_product", "Prokka_inference",  "ICEberg_ID", "ICEberg_Target")
ice <- grepl.sub(gff, "ICE", "feature")
table <- ice[, col]
out <- paste0(output_file, "_ice_gene_table", ".tsv", sep = "")
write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

### Create fields - PHAST
gff$PHAST_ID <- getAttributeField(gff$attributes, "PHAST_ID", ";")
gff$PHAST_Target <- getAttributeField(gff$attributes, "PHAST_Target", ";")

col = c("seqname", "start", "end", "feature", "source", "Prokka_ID", "Prokka_gene", "Prokka_geneFamily", 
        "Prokka_name", "Prokka_product", "Prokka_inference", "PHAST_ID", "PHAST_Target")
prophage <- grepl.sub(gff, "prophage", "feature")
table <- prophage[, col]
out <- paste0(output_file, "_prophage_gene_table", ".tsv", sep = "")
write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Write table
col = c("seqname", "start", "end", "feature", "source", "Prokka_ID", "Prokka_gene", "Prokka_geneFamily", 
        "Prokka_name", "Prokka_product", "Prokka_inference", "VFDB_ID", "VFDB_Target", "Victors_ID", 
        "Victors_Target", "CARD_ARO_Accession", "CARD_Drug_Class", "Resfinder_ID", "Resfinder_Target", 
        "ICEberg_ID", "ICEberg_Target", "PHAST_ID", "PHAST_Target")
table <- gff[, col]
out <- paste0(output_file, "_complete_gene_table", ".tsv", sep = "")
write.table(table, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
