library(rmarkdown)

render("pipeline_report.Rmd", 
       params = list( 
         prefix = "klebsiella_melbourne", 
         fasta = "NB_01/NB01_hybrid_assembly/spades_output/pilon_results/pilon_ERR1023775.fasta", 
         gff = "NB_01/NB01_hybrid_assembly/spades_output/annotation/gffs/final/klebsiella_melbourne_final.gff", 
         gbk = "NB_01/NB01_hybrid_assembly/spades_output/annotation/genbankFile/klebsiella_melbourne.genbank", 
         command = "nextflow run pipeline -c config", 
         config = "NB_01/NB01_hybrid_assembly/spades_output/bacterial_genome_annotation.config", 
         gff_dir = "NB_01/NB01_hybrid_assembly/spades_output/annotation/gffs/final", 
         is_ICEberg = TRUE,
         ice_summary = "NB_01/NB01_hybrid_assembly/spades_output/annotation/summary/klebsiella_melbourne_ices.tsv", 
         ice_gff = "NB_01/NB01_hybrid_assembly/spades_output/annotation/gffs/final/klebsiella_melbourne_ices.gff",
         is_resistance = TRUE,
         resistance_summary = "NB_01/NB01_hybrid_assembly/spades_output/annotation/summary/klebsiella_melbourne_resistance.tsv", 
         resistance_gff = "NB_01/NB01_hybrid_assembly/spades_output/annotation/gffs/final/klebsiella_melbourne_resistance.gff",
         is_CARD = TRUE,
         card_summary = "NB_01/NB01_hybrid_assembly/spades_output/annotation/summary/klebsiella_melbourne_card.tsv",
         is_Resfinder = FALSE,
         is_virulence = TRUE,
         virulence_summary = "NB_01/NB01_hybrid_assembly/spades_output/annotation/summary/klebsiella_melbourne_virulence.tsv",
         virulence_gff = "NB_01/NB01_hybrid_assembly/spades_output/annotation/gffs/final/klebsiella_melbourne_virulence.gff",
         is_VFDB = TRUE,
         is_Victors = TRUE
         ))