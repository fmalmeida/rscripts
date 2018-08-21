#!/usr/bin/Rscript

suppressMessages(library(ggbio))
suppressMessages(library(GenomicRanges))
suppressMessages(library(DataCombine))
suppressMessages(library(ballgown))

f <- function(x, name, value) {
  elementMetadata(x)[[ name ]] <- value
  return(x) }

# Setting parameters
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="gff file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.png",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-b", "--bed"), type = "character", default=NULL,
              help="bed with genome coords [default= %default]", metavar="character")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Getting inputs
gff <- gffRead(opt$input)
gff.gr <- makeGRangesFromDataFrame(gff, keep.extra.columns = TRUE)
genome <- read.delim(opt$bed, header = FALSE,
                     col.names = c("seqname", "start", "end"))
gr <- makeGRangesFromDataFrame(genome)
seqlengths(gr) <- genome$end

# Plot window
window <- 100000L

# ICEberg
ice <- grepl.sub(gff, pattern = "ICE", Var = "feature")
ice.gr <- makeGRangesFromDataFrame(ice, keep.extra.columns = TRUE)
tiles <- tile(gr, width = window)
tiles <- unlist(tiles)
tiles$count <- countOverlaps(tiles, ice.gr, maxgap = 0L, minoverlap = 0L)
tiles.ice.gr <- tiles
tiles.ice.gr$Feature <- "ICEs"

# Virulence
vir <- grepl.sub(gff, pattern = "virulence", Var = "feature")
vir.gr <- makeGRangesFromDataFrame(vir, keep.extra.columns = TRUE)
vir.gr <- f(vir.gr, "Feature", "Virulence")
tiles <- tile(gr, width = window)
tiles <- unlist(tiles)
tiles$count <- countOverlaps(tiles, vir.gr, maxgap = 0L, minoverlap = 0L)
tiles.vir.gr <- tiles
tiles.vir.gr$Feature <- "Virulence"


# Resistance
res <- grepl.sub(gff, pattern = "resistance", Var = "feature")
res.gr <- makeGRangesFromDataFrame(res, keep.extra.columns = TRUE)
res.gr <- f(res.gr, "Feature", "Resistance")
tiles <- tile(gr, width = window)
tiles <- unlist(tiles)
tiles$count <- countOverlaps(tiles, res.gr, maxgap = 0L, minoverlap = 0L)
tiles.res.gr <- tiles
tiles.res.gr$Feature <- "Resistance"

# Phage
phage <- grepl.sub(gff, pattern = "prophage", Var = "feature")
phage.gr <- makeGRangesFromDataFrame(phage, keep.extra.columns = TRUE)
phage.gr <- f(phage.gr, "Feature", "Phage")
tiles <- tile(gr, width = window)
tiles <- unlist(tiles)
tiles$count <- countOverlaps(tiles, phage.gr, maxgap = 0L, minoverlap = 0L)
tiles.phage.gr <- tiles
tiles.phage.gr$Feature <- "Phage"

# Plot
## Getting scales
scale <- round(max(genome$end / as.integer(gsub("L", "", window))))
##

if( length(unique(gff$seqname)) == 1 ) {
  p <- ggbio() +
    circle(tiles.res.gr, geom = "line", aes(color = "Resistance", y=count), radius = 10) +
    circle(tiles.vir.gr, geom = "line", aes(color = "Virulence", y=count), radius = 15) +
    circle(tiles.ice.gr, geom = "line", aes(color = "ICEs", y=count), radius = 20) +
    circle(tiles.phage.gr, geom = "line", aes(color = "Phage", y=count), radius = 25) +
    circle(gr, geom = "ideo", fill = "gray70", radius = 30) +
    circle(gr, geom = "scale", radius = 35, scale.n = scale) +
    circle(gr, geom = "text", aes(label = seqnames), vjust = 0, 
           size = 3.5, radius = 5, angle = 0) +
    theme(legend.position = "bottom")
} else  {
  p <- ggbio() +
    circle(tiles.res.gr, geom = "line", aes(color = "Resistance", y=count), radius = 10) +
    circle(tiles.vir.gr, geom = "line", aes(color = "Virulence", y=count), radius = 15) +
    circle(tiles.ice.gr, geom = "line", aes(color = "ICEs", y=count), radius = 20) +
    circle(tiles.phage.gr, geom = "line", aes(color = "Phage", y=count), radius = 25) +
    circle(gr, geom = "ideo", fill = "gray70", radius = 30) +
    circle(gr, geom = "scale", radius = 35, scale.n = scale) +
    circle(gr, geom = "text", aes(label = seqnames), vjust = 0, size = 3.5, radius = 50, angle = 45) +
    theme(legend.position = "bottom")
}

# Save Image
png <- paste0(opt$out, ".png", collapse = "")
png(png, width = 1280, height = 1080, res = 130)
#png(opt$out)
p # Make plot
dev.off()

svg <- paste0(opt$out, ".svg", collapse = "")
svg(svg, width = 1280, height = 1080, res = 130)
#png(opt$out)
p # Make plot
dev.off()