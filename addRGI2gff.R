#!/usr/bin/Rscript
# Setting Help
'usage: addRGI2gff.R [--input=<file> --gff=<file> --out=<chr> --database=<chr> --type=<chr>]

options:
-g, --gff=<file>      GFF file to add NCBI AMR Annotations into
-i, --input=<file>    RGI tabular output
-o, --out=<chr>       Output file name [default: out.gff]
-t, --type=<chr>      Type of feature. Ex: resistance
-d, --database=<chr>  Name of databased which Blast came from' -> doc

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
#gff <- gffRead(opt$gff)
gff <- gffRead("/work/sample_dataset/annotation/TESTE/resistance/example.gff")
gff$ID <- getAttributeField(as.character(gff$attributes), "ID", ";")

# Load CARD RGI results
rgi_input <- read.delim(
  "/work/sample_dataset/annotation/TESTE/resistance/RGI_annotation/RGI_Example.txt",
  header = TRUE)
## Rename contigs
rgi_input$Contig <- gsub(pattern = "_[0-9]*$", replacement = "", rgi_input$Contig)

if (is.null(rgi_input) == FALSE & dim(rgi_input)[1] != 0) {
  
  rgi_input$Source <- ("CARD-RGI")
  rgi_input$Feature <- ("Resistance")
  rgi_input$Score <- (".")
  rgi_input$Frame <- 0
  rgi_input$Attributes <- 
    paste("card_name=", rgi_input$Best_Hit_ARO, ";rgi_inference=", rgi_input$Model_type,
          ";card_product=", rgi_input$AMR.Gene.Family, ";targeted_drug_class=",
          rgi_input$Drug.Class, ";additional_database=CARD-RGI", sep = "")
  rgi_input$Attributes <- gsub(pattern = " ", replacement = "_",
                               x = rgi_input$Attributes)
  rgi_input$Attributes <- gsub(pattern = "-", replacement = "_",
                               x = rgi_input$Attributes)
         
  rgi_gff <- rgi_input %>% 
    select(Contig, Source, Feature, Start, Stop, 
           Score, Orientation, Frame, Attributes)
  
} else {
  # Load GFF file
  gff <- gffRead(opt$gff)
  # Write output
  write.table(gff, file = opt$out, quote = FALSE, sep = "\t", 
              col.names = FALSE, row.names = FALSE)
}
