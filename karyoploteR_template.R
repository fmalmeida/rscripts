#!/usr/bin/Rscript

################
### Template ###
################

# This document was created to understand how to package works and how to take advantage of it
# to create various karyoplots by simply customizing this file.
#
# Please, do not edit this document directly.
# Copy it to your working directory and work inside your local copy
#
# Felipe Almeida <marques.felipe@aluno.unb.br>


#################
### Tutorials ###
#################

# https://bernatgel.github.io/karyoploter_tutorial/#Tutorial
# https://www.bioconductor.org/packages/release/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html
# https://www.rapidtables.com/web/color/RGB_Color.html

######################
### Load libraries ###
######################

suppressMessages(library(karyoploteR))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))
suppressMessages(library(plyranges))
suppressMessages(library(dplyr))

########################
### Useful functions ###
########################

# Load genome chr from bed (3 cols only)
# It is possible to filter chrs based on min and max sizes
# Default: All chrs
grFromBed <- function(bed, minSize=1, maxSize=NULL) {
  
  if (!is.null(maxSize)) {
    # Open bed
    genome.bed <- read.csv(sep = "\t", header = FALSE,
                           stringsAsFactors = FALSE,
                           file = bed,
                           col.names = c('Chr', 'Start', 'End')) %>%
      dplyr::filter((End - Start) >= minSize & (End - Start) <= maxSize)
    
    # Set BED start to one
    genome.bed$Start <- 1
    
    # Return GRanges
    return(
      with(genome.bed, GRanges(Chr, IRanges(Start, End)))
    )
    
  } else {
    
    # Open bed
    genome.bed <- read.csv(sep = "\t", header = FALSE,
                           stringsAsFactors = FALSE,
                           file = bed,
                           col.names = c('Chr', 'Start', 'End')) %>%
      dplyr::filter((End - Start) >= minSize)
    
    # Set BED start to one
    genome.bed$Start <- 1
    
    # Return GRanges
    return(
      with(genome.bed, GRanges(Chr, IRanges(Start, End)))
    )
    
  }
}

####################
### Loading data ###
####################

# Import GFF as GRanges object
gff.gr <- readGFFAsGRanges("../input/Bexcelsav1.1.gene.gff3")

# Filtering GFF features
## With the function plyranges::filter we can filter a GRanges object based
## on the value of any of its columns
genes.gr <- plyranges::filter(gff.gr, type == "gene")

# Import BED as GRanges
genome.gr <- grFromBed("../input/bnut_chr.bed", minSize = 500000, maxSize = 20e6)

# Importing custom features from a filtered GFF
# We can filter it directly via the GRanges of the complete GFF with the
# function plyranges::filter, or, alternatively, we can filter a GFF and
# create another object directly from the filtered GFF
regions1 <- readGFFAsGRanges("../input/bnut_NLR_annotator_only.gff")

####################
### Drawing Plot ###
####################

## Initiate an instance for the plot
png("teste_type1.png", width = 1200, height = 900)

# Params
## In this area, we use the function getDefaultPlotParams so we have a 
## small list which contains all the necessary parameters used
## when plotting a karyoploteR ideogram. Whith this list we can
## change some parameters and customize the plot a little bit
##
## Read: https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotParams/PlotParams.html
pp <- getDefaultPlotParams(plot.type = 1)
pp$bottommargin <- 500
pp$topmargin <- 500
pp$ideogramheight <- 50
pp$data1outmargin <- 500
pp$data1inmargin <- 200

# Karyotype
## In this area, we use the function plotKaryotype do draw the ideogram
## blocks based on a given GRanges object with the chromosome sizes
##
## kpDataBackground changes the color of the data panel background
## kpAddBaseNumbers is used to draw the position ticks which enable
## the visualization of the genome position (base ticks/markers)
kp <- plotKaryotype(plot.type = 1, main = 'BNUT - Ideogram plot', cex = .9, 
                    genome = genome.gr, plot.params = pp, labels.plotter = NULL)
kpAddChromosomeNames(kp, cex = .5, yoffset = -.01) # The options yoffset and xoffset control the position where the names will appear
kpDataBackground(kp, color = "#F7F7F7")
kpAddBaseNumbers(kp, tick.dist = 2500000, tick.len = 50, cex = .65)

# Gene density
## In this area, we use the function kpPlotDensity to create a histogram
## of the density of genes found in a given genome window size
##
## data.panel = "ideogram" Is used to plot inside the ideograms
kpPlotDensity(kp, genes.gr, window.size = 1e6, 
              data.panel = 1, col="#313D7C", border="#313D7C",
              r0 = 0, r1 = 0.5)
kpAddLabels(kp, labels = "Gene Density", r0=0, r1=0.5, cex = .5)

# Custom regions
## In this are, we use the function kpPlotRegions to draw little blocks
## that represents the genomic regions of all the features inside a given
## GRanges object
##
## kpAddLabels give us a little description of the data
kpPlotRegions(kp, data=regions1, col="#B22222", layer.margin = 0.05, 
              border="#B22222", r0=0.55, r1=1)
kpAddLabels(kp, labels = "Resistance Genes", r0=0.55, r1=1, cex = .5)

# Finish the plot instance
dev.off()
