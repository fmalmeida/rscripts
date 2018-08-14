# Load Library
library(Gviz)
library(DataCombine)
library(ballgown)

# Parameters
gff_file <- "~/Dropbox/prokka_test/out_nextflow/teste_1_merged.gff"
ref_start <- 1
ref_end <- 4700000
ref_chr_name <- 'NC_000913.3'

# Load merged gff as a df
merged_gff <- gffRead(gff_file)
data_g <- GenomicRanges::makeGRangesFromDataFrame(df = merged_gff, seqnames.field = "seqname", start.field = "start", end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
data_track <- AnnotationTrack(data_g, options(ucscChromosomeNames=FALSE), name = "Features", width = 15, showFeatureId = T, min.height=2, feature = "feature")
displayPars(data_track) <- list(col="lightblue", lwd=2, showTitle=TRUE, background.title = "black")

# Set ref track

ref <- GRanges(ref_chr_name, IRanges(ref_start, ref_end))
ref_track <- GenomeAxisTrack(ref, lwd=4, fontsize=20)

# Plot Track
plotTracks(c(ref_track, data_track), main = "Plotting gff features")

# Set especific tracks
resistance_df <- grepl.sub(merged_gff, pattern = "*resistance*", Var = "attributes")
resistance_g <- GenomicRanges::makeGRangesFromDataFrame(df = resistance_df, seqnames.field = "seqname", start.field = "start", end.field = "end", strand.field = "strand", keep.extra.columns = TRUE)
resistance_track <- AnnotationTrack(resistance_g, options(ucscChromosomeNames=FALSE), name = "Resistance", width = 15, showFeatureId = T, min.height=2, feature = "feature")
displayPars(resistance_track) <- list(col="red", lwd=2, showTitle=TRUE, background.title = "black")

plotTracks(c(ref_track, data_track, resistance_track), type = "histogram", main = "Plotting gff features")

# DataTracks

all_features_plus_dttrack <- DataTrack(range = data_g, genome = ref_chr_name, name = "Features in plus strand", strand = "+")

all_features_minus_dttrack <- DataTrack(range = data_g, genome = ref_chr_name, name = "Features in minus strand", strand = "-")

plotTracks(c(ref_track, all_features_plus_dttrack, all_features_minus_dttrack), type = "heatmap")
