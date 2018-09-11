#!/usr/bin/Rscript
# Plot with multiple layers
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

# Setting parameters
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-p", "--pattern"), type = "character", default=NULL,
              help="subset pattern [default= %default]", metavar="character"),
  make_option(c("-f", "--field"), type = "character", default="attributes",
              help="subset pattern [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

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

## Data without NA values
sub <- na.omit(merged_df)

# Count - with NA values
count <- count(merged_df, c("seqname", "ID", "geneFamily", "product"))

id_sb_csv <- paste0(count$seqname, "-", count$geneFamily, 
                    "-", count$product, "-", count$ID,
                    ",", count$freq, sep = "")
write(id_sb_csv, file = "sb.csv")

# plot
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
out <- paste0(opt$out, "_", "with_NA_values", sep = "")
htmlwidgets::saveWidget(sb, opt$out, selfcontained = FALSE)

# Delete temp
file.remove("sb.csv")

# Count - without NA values
count <- count(sub, c("seqname", "ID", "geneFamily", "product"))

id_sb_csv <- paste0(count$seqname, "-", count$geneFamily, 
                    "-", count$product, "-", count$ID,
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
out <- paste0(opt$out, "_", "without_NA_values", sep = "")
htmlwidgets::saveWidget(sb, opt$out, selfcontained = FALSE)

# Delete temp
file.remove("sb.csv")

############################################
## Same plot - Specific for CARD database ##
############################################
gff <- gffRead(opt$input)
card_df <- grepl.sub(gff, pattern = "CARD", Var = "source")

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

## Data without NA values
sub <- na.omit(card_df)

# Count
count <- count(sub, 
               c("seqname", "ID", "GeneFamily", 
                 "DrugClass", "ResistanceMechanism", 
                 "ARO", "Name"))

id_sb_csv <- paste0(count$seqname, "-", count$DrugClass, 
                    "-", count$ResistanceMechanism, "-", 
                    count$GeneFamily, "-", count$Name, 
                    "-", count$ARO, 
                    "-", count$ID, ",", count$freq, sep = "")
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
htmlwidgets::saveWidget(sb, "Complete_CARD_entries.html", 
                        selfcontained = FALSE)

# Delete temp
file.remove("sb.csv")

## Concise CARD plot
id_sb_csv <- paste0(count$seqname, "-", count$DrugClass, 
                    "-", count$ResistanceMechanism, "-", 
                    count$ARO, ",", count$freq, sep = "")
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
htmlwidgets::saveWidget(sb, "Concise_CARD_entries.html", 
                        selfcontained = FALSE)

# Delete temp
file.remove("sb.csv")