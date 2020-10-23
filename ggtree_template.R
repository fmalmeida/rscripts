#############################################################
### Script used for plotting and pruning the phylogenetic ###
### tree of k. variicola genomes produced by parSNP.      ###
###                                                       ###
### Author: Felipe M. Almeida                             ###
#############################################################

###############################
### Load required libraries ###
###############################
library(rentrez)
library(XML)
library(ggplot2)
library(phangorn)
library(ape)
library(ggtree)
library(dplyr)
library(rlist)
library(Polychrome)


########################
### Useful functions ###
########################
source("./r_functions.R") # load function

########################
### Load input files ###
########################
# Input tree
in_tree <- read.tree("./parsnp.tree")
in_tree$tip.label <- shorter_names(in_tree$tip.label)

# Kv15 was used as ref ... we need to change their name
in_tree$tip.label <- gsub(pattern = ".ref", replacement = "", fixed = T,
     x = in_tree$tip.label)

##########################################
### Get metadata -- fetching countries ###
### The for-loop below was used to     ###
### fetch metadata, re-run if you like ###
##########################################

#in_metadata <- data.frame()
#for (entry in in_tree$tip.label) {
#  
#  if (grepl('^Kv', entry)) {
#    country <- "Brazil"
#    isol_source <- "missing"
#    host <- "Homo sapiens"
#    df <- data.frame(entry, country, isol_source, host)
#    colnames(df) <- c("genome", "country", "isolation_source", "host")
#    in_metadata <- rbind(in_metadata, df)
#    
#  } else {
#    in_list <- fetch_metadata(gsub("^\\s+|\\s+$", "", entry))
#    country <- in_list[1]
#    isol_source <- in_list[2]
#    host <- in_list[3]
#    df <- data.frame(entry, country, isol_source, host)
#    colnames(df) <- c("genome", "country", "isolation_source", "host")
#    in_metadata <- rbind(in_metadata, df)
#  }
#  
#}
#write.table(in_metadata, file = "./genome_metadata.txt", quote = FALSE,
#           row.names = FALSE, col.names = FALSE, sep = "\t")

# Input metadata
metadata <- read.csv("./genome_metadata.txt", sep = "\t", header = F, stringsAsFactors = F)
colnames(metadata) <- c("genome", "country", "isolation_source", "host")

######################################################
### parSNP tree roots at the reference genome used ###
### but this is not meaningful to us, thus, let's  ###
### re-root the tree at the midpoint               ###
######################################################
# Re-root tree at midpoint
tree1 <- midpoint(in_tree, node.labels = "support")

# View tree
p <- tree1 %>%
  ggtree(layout = "circular") +
  geom_tippoint() +
  geom_tiplab(align = T)
p


#########################################################
### As we can see, a lot of leaves contain >5 genomes ###
### side-by-side indicating few snps between them.    ###
### This set-up has poor visibility, with extended    ###
### leaves. In order to improve readability we will   ###
### select a few genomes to drop from the tree        ###
#########################################################
# Drop selected genomes from tree and view result
tree1 %>%
  drop.tip(readLines("./ggtree_genomes_for_drop.txt")) %>%
  ggtree(layout='circular') +
  geom_tippoint() +
  geom_tiplab(align = T)

# If good, save it
reduced_tree <- drop.tip(tree1, readLines("./ggtree_genomes_for_drop.txt"))


#####################################################
### Now with a more readable tree, we can add the ###
### countries metadata and plot it in the tree    ###
#####################################################
# Fetch countries
reduced_tree$tip.country <- get_countries(reduced_tree, metadata)

# Transform in data.frame
dd <- do.call(rbind, Map(data.frame, node=1:length(reduced_tree$tip.label), Country=reduced_tree$tip.country))

# Sort countries
dd$Country <- as.character(dd$Country)
dd <- dd[order(dd$Country),]
dd$Country <- as.factor(dd$Country)

# Remove row.names
row.names(dd) <- NULL

# Plot it in the tree
# Showing our genomes and the ones near them
p <- reduced_tree %>%
  ggtree(layout='circular') %<+% dd +
  geom_tippoint(aes(color=Country), size=2, alpha=.9) +
  geom_tiplab(
    size = 2.5,
    linesize = .5,
    align = TRUE,
    aes(subset=label %in% c("Kv15.fna", "Kv35.fna", "Kv57.fna", "Kv97.fna", "Kv104.fna",
                            "GCF_001261855.1", "GCF_004312565.1", "GCF_003255775.1", "GCF_003195165.1"))
  )
p


#################################################
### The colors are not so good. Let's make it ###
### darker so they are better differentiated  ###
#################################################
# Create a custom palette
custom_palette <- 
  createPalette(N = length(unique(reduced_tree$tip.country)),
                seedcolors = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666'))
# Plot with custom palette
p <- reduced_tree %>%
  ggtree(layout='circular') %<+% dd +
  geom_tippoint(aes(color=Country), size=2, alpha=.9) +
  geom_tiplab(
    size = 2.5,
    linesize = .5,
    align = TRUE,
    aes(subset=label %in% c("Kv15.fna", "Kv35.fna", "Kv57.fna", "Kv97.fna", "Kv104.fna",
                            "GCF_001261855.1", "GCF_004312565.1", "GCF_003255775.1", "GCF_003195165.1"))) + 
  scale_color_manual(values = unname(custom_palette))
p

#############################################
### Now it's better. Let's save this plot ###
#############################################
ggsave(filename = "./ggtree_test.png", plot = p, device = "png")
ggsave(filename = "./ggtree_test.svg", plot = p, device = "svg")


################################################################
### Now let's compare the genome content information between ###
### all these genomes ... Using ppanggolin we can query the  ###
### pangenome with a set of genes of interest and discover   ###
### which genomes of the pangenome contain or not the target ###
### gene used as query                                       ###
################################################################
# Load gene presence information as matrix
gene_matrix <- read.table("./gene_presence_absence.Rtab", header = TRUE, stringsAsFactors = F)

# Let's sort this matrix
attach(gene_matrix)
gene_matrix <- gene_matrix[order(Gene),]
detach(gene_matrix)

# We need to transpose this matrix
# in order to plot it
gene_matrix_t <- t(gene_matrix)

# Get colnames (genome names)
colnames(gene_matrix_t) <- gene_matrix_t[1,]
gene_matrix_t <- gene_matrix_t[-1,]

# Small fix on names
rownames(gene_matrix_t) <- shorter_names(rownames(gene_matrix_t)) %>%
  lapply(function(t){
    gsub(x = t, pattern = "gbk", replacement = "fna", fixed = T)
  })

# Change label from number to category
gene_matrix_t[gene_matrix_t==0] <- "absent"
gene_matrix_t[gene_matrix_t==1] <- "present"

# Reload tree
p <- reduced_tree %>%
  ggtree(layout='rectangular', branch.length = "none")

# Plot heatmap
gheatmap(p, gene_matrix_t, colnames_angle = 90, font.size = 3, offset = 1, colnames_offset_y = 15, 
         colnames_position = "top") +
  scale_fill_manual(values=c("#f63737", "#c1eaea"))


#######################################################
### Albeit it shows a nice overview, we can't infer ###
### the distance of the genomes based on the gene   ###
### content. For that, we need to re-calculate the  ###
### tree based on the gene presence matrix as source###
### of distances                                    ###
#######################################################
# We need to transpose this matrix
# in order to plot it
gene_matrix_t <- t(gene_matrix)

# Get colnames (genome names)
colnames(gene_matrix_t) <- gene_matrix_t[1,]
gene_matrix_t <- gene_matrix_t[-1,]

# Small fix on names
rownames(gene_matrix_t) <- shorter_names(rownames(gene_matrix_t)) %>%
  lapply(function(t){
    gsub(x = t, pattern = "gbk", replacement = "fna", fixed = T)
  })

# Calc dist
gene_dist <- dist(gene_matrix_t, method = "manhattan")

# Generate tree
cluster <- hclust(gene_dist)
tree2 <- as.phylo(cluster)

# Check hosts
first_filter <- metadata[metadata$host != "Homo sapiens", ]
second_filter <- first_filter[!grepl("human|urine", first_filter$isolation_source), ]
third_filter <- second_filter[(!grepl("missing", second_filter$isolation_source) | 
                               !grepl("missing", second_filter$host)), ]
final <- third_filter[!grepl("not collected|unknown|respiratory|animal|Chilli|Masala|milk|Chicken", third_filter$isolation_source), ]
from_env <- final$genome

# Reload tree
grp <- list(
  our_samples = c("Kv15.fna", "Kv35.fna", "Kv57.fna", "Kv97.fna", "Kv104.fna"),
  environmental = from_env
)
p <- tree2 %>%
  groupOTU(grp) %>%
  ggtree(layout='rectangular', aes(color=group),
         branch.length = 'none') +
  scale_color_manual(values = c("black", "darkgreen", "firebrick"),
                     name = "Groups",
                     guide = 'none')
p

# Remove some repetitive genomes
sampled <- list.sample(get_taxa_name(p, 478), 50) # Node 478
sampled <- rbind(sampled, list.sample(get_taxa_name(p, 323), 15)) # Node 323
#sampled <- rbind(sampled, list.sample(get_taxa_name(p, 383), 15)) # Node 383
#sampled <- rbind(sampled, list.sample(get_taxa_name(p, 413), 15)) # Node 383

# Drop sampled tips
reduced_tree2 <- tree2 %>%
  drop.tip(sampled) %>%
  groupOTU(grp) %>%
  ggtree(layout='rectangular', aes(color=group),
         branch.length = 'none') +
  scale_color_manual(values = c("black", "darkgreen", "firebrick"),
                     name = "Groups",
                     guide = 'none')

reduced_tree2 + geom_nodepoint()

# Change label from number to category
gene_matrix_t[gene_matrix_t==0] <- "absent"
gene_matrix_t[gene_matrix_t==1] <- "present"

# Plot with heatmap
heat_p <- gheatmap(reduced_tree2, gene_matrix_t, colnames_angle = 90, font.size = 2.5, offset = 1.5, colnames_offset_y = 15, 
         colnames_position = "top") +
  scale_fill_manual(values = c("#f63737", "#c1eaea"),
                    name   = "Virulence genes") +
  theme(
    legend.background=element_rect(fill=NA),
    legend.title=element_text(size=10), # The size should be adjusted with different devout.
    legend.text=element_text(size=8),
    legend.spacing.y = unit(0.25, "cm"))
heat_p

# Save plots
ggsave(filename = "./ggtree_test_virulence_genes.png", plot = heat_p, device = "png")
ggsave(filename = "./ggtree_test_virulence_genes.svg", plot = heat_p, device = "svg")
