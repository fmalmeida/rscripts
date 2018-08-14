# Load Library
library(DataCombine)
library(ballgown)
library(plotly)
library(plyr)

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
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

#WORK
gff <- gffRead(opt$input)

# Create fields
gff$ID <- getAttributeField(gff$attributes, "ID", ";")
gff$gene <- getAttributeField(gff$attributes, "gene", ";")
gff$name <- getAttributeField(gff$attributes, "Name", ";")
gff$product <- getAttributeField(gff$attributes, "product", ";")

# Write table
table <- gff[,-c(1:9)]

# Filter
sub_df <- grepl.sub(gff, pattern = "*resistance*", Var = "attributes")
card_df <- grepl.sub(gff, pattern = "*resistance*", Var = c("feature"))

merdeg_df <- merge.data.frame(sub_df, card_df, all = TRUE)

merdeg_df$ID <- getAttributeField(merdeg_df$attributes, "ID", ";")
merdeg_df$gene <- getAttributeField(merdeg_df$attributes, "gene", ";")
merdeg_df$name <- getAttributeField(merdeg_df$attributes, "Name", ";")
merdeg_df$product <- getAttributeField(merdeg_df$attributes, "product", ";")
merdeg_df$product <- sub("Multiple_antibiotic_resistance_protein", "Multidrug_resistance_protein", 
                         merdeg_df$product)
merdeg_df$product <- sub("_[^_]+$", "", merdeg_df$product)
merdeg_df$product <- gsub("-", "_", merdeg_df$product)

#Count
count <- count(merdeg_df, c("seqname", "product"))


# Plot
p <- plot_ly(count, labels = ~product, values = ~freq, type = "pie", 
             textposition = 'outside',
             textinfo = 'label+value') %>%
  layout(title = 'Count of resistance products', showlegend = FALSE,
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

# Save
htmlwidgets::saveWidget(as_widget(p), file = opt$out)
