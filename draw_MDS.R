#!/usr/bin/Rscript
# A supporting script for drawing the multidimensional scalling analyses!
# To be used with fmalmeida/BacterialMDS
# Author: Felipe Marques de Almeida <almeidafmarques@gmail.com>

#################################################
### Setting Help message and input parameters ###
#################################################
'usage: draw_MDS.R [--input=<file>]
options:
-i, --input=<file>    Quadratic distance matrix produced with Dashing software
-m, --meta=<file>     Path to file containing metadata information of genomes
-q, -qmeta=<file>     Path to file containing metadata information of query genomes' -> doc

# Parse parameters
suppressMessages(library(docopt))
opt <- docopt(doc)

if (is.null(opt$input)){
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

#####################################
### Loading of necessary packages ###
#####################################
library(MASS)
library(dplyr)
library(vegan)
library(ggplot2)
library(stringr)
library(plotly)
library(tidyverse)
library(htmlwidgets)

#############################
######### Functions #########
#############################

## Create function for loading matrices
load_dist_matrix <- function(matrix) {
  
  # Load matrix
  matrix_df <- read.csv(matrix, sep = "\t", header = FALSE)
  
  # Set rownames
  rownames(matrix_df) <- matrix_df[,1]
  
  # Remove unwanted col
  matrix_df <- dplyr::select(matrix_df, -V1)
  
  # Set as dist
  dist <- as.dist(as(matrix_df, "matrix"))
  
  return(dist)
  
}

## Create function to get accessions
get_acc <- function(df) {
  splits <- strsplit(sapply(strsplit(as.character(df$X.Names), "\\/"), `[`, 2), "_", fixed = TRUE)
  
  return(
    lapply(splits,FUN=function(x){
      paste(x[1],x[2],sep="_")
    })
  )
}

## Subset a string and get the last two values
substr_latest <- function(x, sep) {
  
  vals <- str_split(x, pattern = sep)
  final <- paste0(tail(vals[[1]], 2), collapse = "~~~")
  
  return(final)
}

## Fix Group Names for shapes and colors
fix_groups <- function(x) {
  if (x[5] == "From database") {
    group <- x[6]
  } else {
    group <- x[5]
  }
  
  return(group[[1]])
}

#################################
######### Load the data #########
#################################

## Loading of distance matrices
dist.mat <- load_dist_matrix("/Volumes/falmeida1TB/Git_Repos/BacterialMDS/output/dashing_distance/dashing_distance_matrix.txt")

###################################
######### Do MDS analysis #########
###################################

## Perform MDS using the Classical Metric
dist.mds <- cmdscale(dist.mat, eig=TRUE)

## Get the exact points of the MDS produced
dist.points <- as.data.frame(dist.mds$points)

## Include names of samples/genomes
dist.points$names <- rownames(dist.points)
dist.points$names <- lapply(dist.points$names, function(x) {
  substr_latest(x, "/")
})

## Getting Identification from file names
splits <- strsplit(sapply(strsplit(as.character(dist.points$names), "~~~"), tail, 1), "_", fixed = TRUE)
dist.points$Accession <- 
  lapply(splits,FUN=function(x){
    if (is.na(x[2])) {
      paste(x[1])
    } else {
      paste(x[1], x[2],sep="_")
    }
  })

## Copy
dist.points$Identification <- dist.points$Accession

## Parse Species / Group
dist.points$Group <- sapply(strsplit(as.character(dist.points$names), "~~~"), head, 1)

## Let only the names of our samples in the data.frame
## To make subsetting and plotting easier
dist.points[!dist.points$Group %in% 
            c("queries"), "Identification"] <- "From database"

## Fix for plots
rownames(dist.points) <- NULL
dist.points$Identification <- apply(dist.points, 1, fix_groups)


#################################
######### Load metadata #########
#################################

## Loading metadata (only of the closest genomes)
metadata <- read.csv("/Volumes/falmeida1TB/Git_Repos/BacterialMDS/output/refseq/Klebsiella_quasipneumoniae_genomes_metadata.txt",
                        sep = "\t",
                        header = FALSE,
                        col.names = c("Assembly.acc", "Biosample.acc", "Owned.by", "Date", "Geo.loc", "Isolation.source", "Host"),
                        stringsAsFactors = FALSE)

query_metadata <- read.csv("/Volumes/falmeida1TB/Git_Repos/BacterialMDS/query_metadata.csv",
                     header = FALSE,
                     col.names = c("Assembly.acc", "Biosample.acc", "Owned.by", "Date", "Geo.loc", "Isolation.source", "Host"),
                     stringsAsFactors = FALSE)

metadata <- rbind(metadata, query_metadata)
metadata$Date <- str_replace_all(string = metadata$Date, "-|_", "/")
splits <- strsplit(sapply(strsplit(as.character(metadata$Assembly.acc), "~~~"), tail, 1), "_", fixed = TRUE)
metadata$Assembly.acc <- 
  lapply(splits,FUN=function(x){
    if (is.na(x[2])) {
      paste(x[1])
    } else {
      paste(x[1], x[2],sep="_")
    }
})

## Merging with metadata
dist.points <- merge(dist.points, metadata,
                   by.x = "Accession", by.y = "Assembly.acc",
                   sort = FALSE, all.x = TRUE, all.y = TRUE)


########################################
######### Begin of plot config #########
########################################

## Dataframe with only our samples (For custom colouring)
queries <- dist.points[dist.points$Group %in% c("queries"), ]

## Setting up the plot
full_plot <- 
  ggplot(dist.points, aes(x=V1, y=V2, shape=Group, color=factor(unlist(Identification)),
                        text = paste(
                          "Owned by: ", Owned.by, "\n",
                          "Accession: ", Accession, "\n",
                          "Biosample: ", Biosample.acc, "\n",
                          "Date: ", Date, "\n",
                          "Isolation source: ", Isolation.source, "\n",
                          "Host: ", Host, "\n",
                          "From: ", Geo.loc, "\n",
                          sep = ""
                        ))) + 
  xlim(min(dist.points$V1), max(dist.points$V1)) +
  ylim(min(dist.points$V2), max(dist.points$V2)) +
  scale_color_brewer(palette="Dark2") +
  geom_point(alpha=0.45) +
  geom_point(data=queries, aes(x=V1, y=V2, shape=Group, color=factor(unlist(Identification)))) +
  theme_classic() + 
  theme(legend.position = "right") +
  labs(y= "", x = "", color = "Colors", shape = "Shapes") +
  ggtitle("Multidimensional scalling of genetic distances")

## Checking
full_plot


##########################################
######### Check interactive plot #########
##########################################
interactive_plot <- ggplotly(full_plot, tooltip = "text")
withr::with_dir(".", htmlwidgets::saveWidget(as_widget(interactive_plot), "interactive_mds_of_genetic_distances.html", selfcontained = TRUE))


#################################
######### Plots outputs #########
#################################

# Saving plots
## SVG
ggsave(filename = "./plots/mds_of_genetic_distances.svg", plot = full_plot, device = "svg", dpi = 1080, width=7, height=7)

## PNG
ggsave(filename = "./plots/mds_of_genetic_distances.png", plot = full_plot, device = "png", dpi = 1080, width=7, height=7)
