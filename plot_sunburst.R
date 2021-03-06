#!/usr/bin/Rscript
# Setting help
'usage: plot_sunburst.R [--input=<file> --out=<chr> --pattern=<chr> --field=<chr>]

options:
  -i, --input=<file>    GFF from which you want to plot
  -o, --out=<chr>       Output prefix file name [default: out]
  -p, --pattern=<chr>   Pattern that will be used to subset GFF [default: CARD]
  -f, --field=<chr>     GFF field to search for pattern [default: source]' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)
if (is.null(opt$input)){
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

# Load libraries
suppressMessages(library(sunburstR))
suppressMessages(library(plyr))
suppressMessages(library(DataCombine))
suppressMessages(library(ballgown))

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

if (opt$pattern == "CARD" && opt$field == "source") {
############################################
## Same plot - Specific for CARD database ##
############################################
gff <- gffRead(opt$input)
card_df <- grepl.sub(gff, pattern = "CARD", Var = "source")

if (is.data.frame(card_df) && nrow(card_df)!=0) {

card_df$ARO <- getAttributeField(card_df$attributes, "ARO", ";")

card_df$DrugClass <- getAttributeField(card_df$attributes, 
                                       "Drug_Class", ";")
card_df$DrugClass <- gsub("-", "_", card_df$DrugClass)
card_df$DrugClass <- gsub(",", "_", card_df$DrugClass)

card_df$GeneFamily <- getAttributeField(card_df$attributes, 
                                        "Gene_Family", ";")
card_df$GeneFamily <- gsub("-", "_", card_df$GeneFamily)
card_df$GeneFamily <- gsub(",", "_", card_df$GeneFamily)

card_df$ResistanceMechanism <- getAttributeField(card_df$attributes, 
                                                 "Resistance_Mechanism", ";")
card_df$ResistanceMechanism <- gsub("-", "_", card_df$ResistanceMechanism)
card_df$ResistanceMechanism <- gsub(",", "_", card_df$ResistanceMechanism)

card_df$Name <- getAttributeField(card_df$attributes, 
                                  "DB_Name", ";")
card_df$Name <- gsub("-", "_", card_df$Name)
card_df$Name <- gsub(",", "_", card_df$Name)
card_df$ID <- getAttributeField(card_df$attributes, "ID", ";")

## Removing NA values
sub <- na.omit(card_df)

# Count
count <- count(sub, 
               c("seqname", "DrugClass", "ResistanceMechanism", "GeneFamily", "ID", "ARO"))

## Concise CARD plot
id_sb_csv <- paste0(count$seqname, "-", count$DrugClass, 
                    "-", count$ResistanceMechanism, "-", 
                    count$GeneFamily, "-", count$ID, ",", count$freq, sep = "")
write(id_sb_csv, file = "sb.csv")
    
#plot
df_sb <- read.delim("sb.csv", header = FALSE, sep = ",", 
                    dec = NULL, quote = NULL)
sb <- sunburst(
df_sb,
count = TRUE, # add count just for demonstration
legend = list(w=200, h=50, r=10), # make extra room for our legend
legendOrder = c(unique(vapply(strsplit(as.character(df_sb[,1]),"-"), 
                              `[`, 1, FUN.VALUE=character(1))), 
                unique(vapply(strsplit(as.character(df_sb[,1]),"-"), 
                              `[`, 2, FUN.VALUE=character(1))))
)
    
# Save
htmlwidgets::saveWidget(sb, paste0(opt$out, "_CARD_compliant.html", sep = ""), 
                            selfcontained = FALSE)
    
# Delete temp
file.remove("sb.csv") } else {
      
}} else {
# WORK
gff <- gffRead(opt$input)

# Create fields
gff$attributes <- gsub(",", ";", gff$attributes)
gff$ID <- getAttributeField(gff$attributes, "ID", ";")
gff$gene <- getAttributeField(gff$attributes, "gene", ";")
gff$geneFamily <- substr(gff$gene, 1, 3)
gff$product <- getAttributeField(gff$attributes, "product", ";")

# Filter
merged_df <- grepl.sub(gff, pattern = opt$pattern, Var = opt$field)
if (is.data.frame(merged_df) && nrow(merged_df)!=0) {

## Altering some tag problems
merged_df$ID <- getAttributeField(merged_df$attributes, "ID", ";")
merged_df$gene <- getAttributeField(merged_df$attributes, "gene", ";")
merged_df$name <- getAttributeField(merged_df$attributes, "Name", ";")
merged_df$product <- getAttributeField(merged_df$attributes, "product", ";")
merged_df$product <- sub("Multiple_antibiotic_resistance_protein", "Multidrug_resistance_protein", 
                         merged_df$product)
merged_df$product <- sub("_[^_]+$", "", merged_df$product)
merged_df$product <- gsub("-", "_", merged_df$product)
merged_df$product <- gsub(",", "_", merged_df$product)

## Removing NA values
sub <- na.omit(merged_df)

# Count feature repetition
count <- count(sub, c("seqname", "product", "name", "ID"))

id_sb_csv <- paste0(count$seqname, "-", count$product, 
                    "-", count$name, "-", count$ID,
                    ",", count$freq, sep = "")
write(id_sb_csv, file = "sb.csv")

#plot
df_sb <- read.delim("sb.csv", header = FALSE, sep = ",", 
                    dec = NULL, quote = NULL)
sb <- sunburst(
  df_sb,
  count = TRUE, # add count just for demonstration
  legend = list(w=200, h=50, r=10), # make extra room for our legend
  legendOrder = c(unique(vapply(strsplit(as.character(df_sb[,1]),"-"), 
                                `[`, 1, FUN.VALUE=character(1))), 
                  unique(vapply(strsplit(as.character(df_sb[,1]),"-"), 
                                `[`, 2, FUN.VALUE=character(1))))
)

# Save
htmlwidgets::saveWidget(sb, paste0(opt$out, ".html", sep = ""), 
                        selfcontained = FALSE)

# Delete temp
file.remove("sb.csv") } else {
  
}}