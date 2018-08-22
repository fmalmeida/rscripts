#!/usr/bin/Rscript
suppressMessages(library(Gviz))
suppressMessages(library(GenomicRanges))
suppressMessages(library(ballgown))
suppressMessages(library(DataCombine))

# Teste para plotar features com outros plots além do circular

gff <- gffRead("~/Documents/testes/output_teste/gffs/merged/ncbi_final.gff")
gff.gr <- makeGRangesFromDataFrame(gff, keep.extra.columns = FALSE)

##Resistance Features
gff.gr <- makeGRangesFromDataFrame(gff, keep.extra.columns = FALSE)
resistance <- grepl.sub(gff, "resistance", "feature")
res.gr <- makeGRangesFromDataFrame(resistance)
strand(res.gr) <- "*"
gff.gr$resistance <-countOverlaps(gff.gr, res.gr)
subR.gr <- subsetByOverlaps(gff.gr, res.gr)
strand(subR.gr) <- "*"

##Virulence Features
gff.gr <- makeGRangesFromDataFrame(gff, keep.extra.columns = FALSE)
virulence <- grepl.sub(gff, "virulence", "feature")
vir.gr <- makeGRangesFromDataFrame(virulence)
strand(vir.gr) <- "*"
gff.gr$virulence <-countOverlaps(gff.gr, vir.gr)
subV.gr <- subsetByOverlaps(gff.gr, vir.gr)
strand(subV.gr) <- "*"

##(ICEs) Features
gff.gr <- makeGRangesFromDataFrame(gff, keep.extra.columns = FALSE)
ices <- grepl.sub(gff, "ICE", "feature")
ices.gr <- makeGRangesFromDataFrame(ices)
strand(ices.gr) <- "*"
gff.gr$ices <-countOverlaps(gff.gr, ices.gr)
subI.gr <- subsetByOverlaps(gff.gr, ices.gr)
strand(subI.gr) <- "*"

##Plotting
rtrack <- DataTrack(subR.gr, options(ucscChromosomeNames=FALSE),
                    name = "Resistance", type=c("heatmap"))
vtrack <- DataTrack(subV.gr, options(ucscChromosomeNames=FALSE),
                    name = "Virulence", type=c("heatmap"))
itrack <- DataTrack(subI.gr, options(ucscChromosomeNames=FALSE),
                    name = "ICEbergs", type=c("heatmap"))
gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack, rtrack, vtrack, itrack), 
           main = "Genomic Features")

# Save Image
png <- paste0("~/Downloads/teste_without_window_and_not_separeted_grs", 
              ".png", collapse = "")
png(png, width = 1280, height = 1080, res = 130)
plotTracks(list(gtrack, rtrack, vtrack, itrack), 
           main = "Genomic Features")
dev.off()

svg <- paste0("~/Downloads/teste_without_window_and_not_separeted_grs", 
              ".svg", collapse = "")
svg(svg, width = 14, height = 7)
plotTracks(list(gtrack, rtrack, vtrack, itrack), 
           main = "Genomic Features")
dev.off()