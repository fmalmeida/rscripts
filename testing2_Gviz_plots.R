#!/usr/bin/Rscript
suppressMessages(library(Gviz))
suppressMessages(library(GenomicRanges))
suppressMessages(library(ballgown))
suppressMessages(library(DataCombine))
# Teste para plotar features com outros plots além do circular

gff <- gffRead("~/Documents/testes/output_teste/gffs/merged/ncbi_final.gff")
gff.gr <- makeGRangesFromDataFrame(gff, keep.extra.columns = FALSE)

# Plot window
window <- 100000L

#BED
genome <- data.frame("NZ_CP027060.1", 1, 4700000)
colnames(genome) <- c("seqname", "start", "end")
gr <- makeGRangesFromDataFrame(genome)
seqlengths(gr) <- genome$end

tiles <- tile(gr, width = window)
tiles <- unlist(tiles)

# ICEberg
ice <- grepl.sub(gff, pattern = "ICE", Var = "feature")
ice.gr <- makeGRangesFromDataFrame(ice, keep.extra.columns = FALSE)
tiles$countIce <- countOverlaps(tiles, ice.gr, maxgap = 0L, minoverlap = 0L)

# Virulence
vir <- grepl.sub(gff, pattern = "virulence", Var = "feature")
vir.gr <- makeGRangesFromDataFrame(vir, keep.extra.columns = FALSE)
tiles$countVir <- countOverlaps(tiles, vir.gr, maxgap = 0L, minoverlap = 0L)


# Resistance
res <- grepl.sub(gff, pattern = "resistance", Var = "feature")
res.gr <- makeGRangesFromDataFrame(res, keep.extra.columns = FALSE)
tiles$countRes <- countOverlaps(tiles, res.gr, maxgap = 0L, minoverlap = 0L)

# Phage
phage <- grepl.sub(gff, pattern = "prophage", Var = "feature")
phage.gr <- makeGRangesFromDataFrame(phage, keep.extra.columns = FALSE)
tiles$countPhage <- countOverlaps(tiles, phage.gr, maxgap = 0L, minoverlap = 0L)


##Plotting

### Generating commentary

sub <- paste0("This heatmap was rendered using a Genome window size of ", 
              window, " base pairs\n", 
              "Each block represent such interval.", collapse = "" )

cTrack <- CustomTrack(plottingFunction = 
                        function(GdObject, prepare=FALSE) 
                        {if(!prepare) grid.text(sub); 
                          return(invisible(GdObject))}, 
                      rot.title = 0, background.title = "transparent", 
                      name = "Comments")

### Really plotting
dtrack <- DataTrack(tiles, options(ucscChromosomeNames=FALSE), 
                    groups = c("countIce", "countVir", "countRes", "countPhage"),
                    name = "Features", type=c("heatmap"), background.title = "black")

gtrack <- GenomeAxisTrack(add53 = TRUE, add35 = TRUE, 
                          col = "black", fontcolor = "black")
plotTracks(list(gtrack, dtrack, cTrack), 
           main = "Genomic Features", box.legend = TRUE)

# Save Image
png <- paste0("~/Downloads/teste_with_window_and_not_separeted_grs", 
              ".png", collapse = "")
png(png, width = 1280, height = 1080, res = 130)
plotTracks(list(gtrack, dtrack, cTrack), 
           main = "Genomic Features", box.legend = TRUE)
dev.off()

svg <- paste0("~/Downloads/teste_with_window_and_not_separeted_grs", 
              ".svg", collapse = "")
svg(svg, width = 14, height = 7)
plotTracks(list(gtrack, dtrack, cTrack), 
           main = "Genomic Features", box.legend = TRUE)
dev.off()

#####################
## Yet another try ##
#####################

tiles <- tile(gr, width = window)
tiles <- unlist(tiles)

# ICEberg
ice <- grepl.sub(gff, pattern = "ICE", Var = "feature")
ice.gr <- makeGRangesFromDataFrame(ice, keep.extra.columns = FALSE)
tiles.ice.gr <- subsetByOverlaps(tiles, ice.gr, maxgap = 0L, minoverlap = 0L)
tiles.ice.gr$count <- countOverlaps(tiles.ice.gr, ice.gr, maxgap = 0L, minoverlap = 0L)

# Virulence
vir <- grepl.sub(gff, pattern = "virulence", Var = "feature")
vir.gr <- makeGRangesFromDataFrame(vir, keep.extra.columns = FALSE)
tiles.vir.gr <- subsetByOverlaps(tiles, vir.gr, maxgap = 0L, minoverlap = 0L)
tiles.vir.gr$count <- countOverlaps(tiles.vir.gr, vir.gr, maxgap = 0L, minoverlap = 0L)

# Resistance
res <- grepl.sub(gff, pattern = "resistance", Var = "feature")
res.gr <- makeGRangesFromDataFrame(res, keep.extra.columns = FALSE)
tiles.res.gr <- subsetByOverlaps(tiles, res.gr, maxgap = 0L, minoverlap = 0L)
tiles.res.gr$count <- countOverlaps(tiles.res.gr, res.gr, maxgap = 0L, minoverlap = 0L)

# Phage
phage <- grepl.sub(gff, pattern = "prophage", Var = "feature")
phage.gr <- makeGRangesFromDataFrame(phage, keep.extra.columns = FALSE)
tiles.phage.gr <- subsetByOverlaps(tiles, phage.gr, maxgap = 0L, minoverlap = 0L)
tiles.phage.gr$count <- countOverlaps(tiles.phage.gr, phage.gr, maxgap = 0L, minoverlap = 0L)

### Plot
rtrack <- DataTrack(tiles.res.gr, options(ucscChromosomeNames=FALSE),
                    name = "Resistance", type=c("heatmap"), background.title = "black")
vtrack <- DataTrack(tiles.vir.gr, options(ucscChromosomeNames=FALSE),
                    name = "Virulence", type=c("heatmap"), background.title = "black")
itrack <- DataTrack(tiles.ice.gr, options(ucscChromosomeNames=FALSE),
                    name = "ICEbergs", type=c("heatmap"), background.title = "black")
ptrack <- DataTrack(tiles.phage.gr, options(ucscChromosomeNames=FALSE),
                    name = "Prophages", type=c("heatmap"), background.title = "black")

gtrack <- GenomeAxisTrack(add53 = TRUE, add35 = TRUE, 
                          col = "black", fontcolor = "black")
plotTracks(list(gtrack, rtrack, vtrack, itrack, ptrack), main = "Genomic Features", 
           sizes = NULL, box.legend = TRUE)


# Save Image
png <- paste0("~/Downloads/teste_with_window_and_separeted_grs", 
              ".png", collapse = "")
png(png, width = 1280, height = 1080, res = 130)
#png(opt$out)
plotTracks(list(gtrack, rtrack, vtrack, itrack, ptrack, cTrack), main = "Genomic Features", 
           sizes = NULL, box.legend = TRUE)
dev.off()

svg <- paste0("~/Downloads/teste_with_window_and_separeted_grs", 
              ".svg", collapse = "")
svg(svg, width = 14, height = 7)
plotTracks(list(gtrack, rtrack, vtrack, itrack, ptrack, cTrack), main = "Genomic Features", 
           sizes = NULL, box.legend = TRUE)
dev.off()

## With same limit

### Plot
rtrack <- DataTrack(tiles.res.gr, options(ucscChromosomeNames=FALSE),
                    name = "Resistance", type=c("heatmap"), 
                    background.title = "black", ylim = c(0,40))
vtrack <- DataTrack(tiles.vir.gr, options(ucscChromosomeNames=FALSE),
                    name = "Virulence", type=c("heatmap"), 
                    background.title = "black", ylim = c(0,40))
itrack <- DataTrack(tiles.ice.gr, options(ucscChromosomeNames=FALSE),
                    name = "ICEbergs", type=c("heatmap"), 
                    background.title = "black", ylim = c(0,40))
ptrack <- DataTrack(tiles.phage.gr, options(ucscChromosomeNames=FALSE),
                    name = "Prophages", type=c("heatmap"), 
                    background.title = "black", ylim = c(0,40))

gtrack <- GenomeAxisTrack(add53 = TRUE, add35 = TRUE, 
                          col = "black", fontcolor = "black")
png <- paste0("~/Downloads/teste_with_window_same_limit_and_separeted_grs", 
              ".png", collapse = "")
png(png, width = 1280, height = 1080, res = 130)
plotTracks(list(gtrack, rtrack, vtrack, itrack, ptrack, cTrack), main = "Genomic Features", 
           sizes = NULL, box.legend = TRUE)
dev.off()

svg <- paste0("~/Downloads/teste_with_window_same_limit_and_separeted_grs", 
              ".svg", collapse = "")
svg(svg, width = 14, height = 7)
plotTracks(list(gtrack, rtrack, vtrack, itrack, ptrack, cTrack), main = "Genomic Features", 
           sizes = NULL, box.legend = TRUE)
dev.off()