#!/usr/bin/Rscript
# Setting Help
'usage: addRGI2gff.R [--input=<file> --gff=<file> --out=<chr>]

options:
-g, --gff=<file>      GFF file to add NCBI AMR Annotations into
-i, --input=<file>    RGI tabular output
-o, --out=<chr>       Output file name [default: out.gff]' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)

if (is.null(opt$gff)){
  stop("At least one argument must be supplied (gff file)\n", call.=FALSE)
}

if (is.null(opt$input)){
  stop("At least one argument must be supplied (AMRFinder output file)\n", call.=FALSE)
}

# Load libraries
suppressMessages(library(ballgown))
suppressMessages(library(DataCombine))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))

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

# Create a column in gff with ids
gff$ID <- getAttributeField(gff$attributes, "ID", ";")

# Load CARD RGI results
rgi_input <- read.delim(opt$input, header = TRUE)

if (is.null(rgi_input) == FALSE & dim(rgi_input)[1] != 0) {
  
  # Drop unused ID column
  rgi_input <- rgi_input %>%
    select(-ID)
  
  # Fix ORF_ID
  rgi_input <- rgi_input %>% mutate(ORF_ID = str_replace(ORF_ID, "\\s", "|")) %>%
    separate(ORF_ID, into = c("Protein_ID", "Protein_Product"), sep = "\\|")
  
  # Fix Drug Classes and Resistance mechanism
  rgi_input$Drug.Class <- gsub(pattern = ";\\s", replacement = "/",
                               x = rgi_input$Drug.Class)
  rgi_input$Resistance.Mechanism <- gsub(pattern = ";\\s", replacement = "/",
                               x = rgi_input$Resistance.Mechanism)
  
  # Remove low identity hits
  rgi_input <- rgi_input[rgi_input$Best_Identities >= 85, ]
  
  # Create Attributes
  rgi_input$NEW_attributes <-
    paste("CARD_ID=", rgi_input$ARO, ";CARD_Target=", rgi_input$Best_Hit_ARO, ";RGI_inference=", rgi_input$Model_type,
          ";CARD_Resistance_Mechanism=", rgi_input$Resistance.Mechanism, ";CARD_drug_class=",
          rgi_input$Drug.Class, ";Additional_database=CARD", sep = "")
  
  # Get ids
  ids <- rgi_input$Protein_ID
  
  # Subset based on gene IDs
  sub <- gff %>% filter(ID %in% ids) %>% select(seqname, source, feature, start, end, score, strand, frame, attributes, ID)
  not <- gff %>% filter(ID %ni% ids) %>% select(seqname, source, feature, start, end, score, strand, frame, attributes)
  
  # Change fields values
  ## source
  s <- sub$source
  sn <- "CARD"
  snew <- paste(s, sn, sep = ",")
  sub$source <- snew
  
  ## feature
  f <- sub$feature
  fn <- "Resistance"
  fnew <- paste(f, fn, sep = ",")
  sub$feature <- fnew
  
  ## attributes
  sub <- merge.data.frame(sub, rgi_input, by.x = "ID",
                          by.y = "Protein_ID", all = TRUE)
  sub <- unite(sub, "attributes", c("attributes", "NEW_attributes"), sep = ";") %>%
    select(seqname, source, feature, start, end, score, strand, frame, attributes)
  
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
