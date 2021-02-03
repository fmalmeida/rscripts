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
in_tree$tip.label <- gsub(pattern = ".fna", replacement = "", fixed = T,
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
    aes(subset=label %in% c("Kv15", "Kv35", "Kv57", "Kv97", "Kv104",
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
    aes(subset=label %in% c("Kv15", "Kv35", "Kv57", "Kv97", "Kv104",
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
gene_matrix <- gene_matrix[order(Gene_Class, Gene),]
detach(gene_matrix)

# We need to transpose this matrix
# in order to plot it
gene_matrix_t <- gene_matrix %>% select(-Gene_Class) %>% t()

# Get colnames (genome names)
colnames(gene_matrix_t) <- gene_matrix_t[1,]
gene_matrix_t <- gene_matrix_t[-1,]

# Small fix on names
rownames(gene_matrix_t) <- shorter_names(rownames(gene_matrix_t)) %>%
  lapply(function(t){
    gsub(x = t, pattern = "gbk", replacement = "fna", fixed = T)
    gsub(x = t, pattern = ".gbk", replacement = "", fixed = T)
  })

# Remove genomes not used before
gene_matrix_t <- gene_matrix_t[row.names(gene_matrix_t) %in% as.list(reduced_tree$tip.label), ]

# Calc dist
gene_dist <- dist(gene_matrix_t, method = "euclidean")

# Generate tree
cluster <- hclust(gene_dist, method = "single")
tree2 <- as.phylo(cluster)

# Reload tree
grp <- list(
  our_samples = c("Kv15", "Kv35", "Kv57", "Kv97", "Kv104")
)

# Drop tips from the bigger leaf (with lots of repetition)
reduced_tree2 <- drop.tip(tree2, readLines("./ggtree_genomes_to_drop_for_heatmap.txt")) %>%
  groupOTU(grp)

# Check results
p <- reduced_tree2 %>%
  ggtree(layout='rectangular', aes(color=group), branch.length = 'none') +
  scale_color_manual(values = c("black", "firebrick"),
                     name = "Groups",
                     guide = 'none') +
  geom_tiplab(
    size = 2.5,
    linesize = .5,
    align = TRUE,
    aes(subset=label %in% c("Kv15", "Kv35", "Kv57", "Kv97", "Kv104")))
p

# Change label from number to category
gene_matrix_t[gene_matrix_t==0] <- "absent"
gene_matrix_t[gene_matrix_t==1] <- "present"

# Change colors by gene class
classes <- gene_matrix %>% select(Gene_Class) %>% unique()
for (i in as.list(classes$Gene_Class)) {
  genes <- gene_matrix %>% filter(Gene_Class == i) %>% select(Gene) %>% unique()
  
  for (u in as.list(genes$Gene)) {
    gene_matrix_t[gene_matrix_t[, u] == "present", u] <- 
      gsub("present", i, gene_matrix_t[gene_matrix_t[, u] == "present", u])
  }

}

# Plot with heatmap
heat_p <- gheatmap(p, gene_matrix_t, colnames_angle = 45, font.size = 2.5, offset = 10, colnames_offset_y = 0.5, 
         colnames_position = "top", color = "black", hjust = 0, width = 2) +
  scale_fill_manual(values=c('white', 'lightgray', '#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69'),
                    name   = "Virulence genes") +
  theme(
    legend.background=element_rect(fill=NA),
    legend.title=element_text(size=10), # The size should be adjusted with different devout.
    legend.text=element_text(size=8),
    legend.spacing.y = unit(0.25, "cm"))
heat_p

# Save plots
ggsave(filename = "./ggtree_virulence_genes_colored.png", plot = heat_p, device = "png")
ggsave(filename = "./ggtree_virulence_genes_colored.svg", plot = heat_p, device = "svg")

###########################################
### Plot virulence genes with countries ###
###########################################
# Load gene presence information as matrix
gene_matrix <- read.table("./gene_presence_absence.Rtab", header = TRUE, stringsAsFactors = F)

# Let's sort this matrix
attach(gene_matrix)
gene_matrix <- gene_matrix[order(Gene_Class, Gene),]
detach(gene_matrix)

# We need to transpose this matrix
# in order to plot it
gene_matrix_t <- gene_matrix %>% select(-Gene_Class) %>% t()

# Get colnames (genome names)
colnames(gene_matrix_t) <- gene_matrix_t[1,]
gene_matrix_t <- gene_matrix_t[-1,]

# Small fix on names
rownames(gene_matrix_t) <- shorter_names(rownames(gene_matrix_t)) %>%
  lapply(function(t){
    gsub(x = t, pattern = "gbk", replacement = "fna", fixed = T)
    gsub(x = t, pattern = ".gbk", replacement = "", fixed = T)
  })

# Remove genomes not used before
gene_matrix_t <- gene_matrix_t[row.names(gene_matrix_t) %in% as.list(reduced_tree$tip.label), ]

# Calc dist
gene_dist <- dist(gene_matrix_t, method = "euclidean")

# Generate tree
cluster <- hclust(gene_dist, method = "single")
tree3 <- as.phylo(cluster)

# Reload tree
grp <- list(
  our_samples = c("Kv15", "Kv35", "Kv57", "Kv97", "Kv104")
)

# Drop tips from the bigger leaf (with lots of repetition)
reduced_tree3 <- drop.tip(tree3, readLines("./ggtree_genomes_to_drop_for_heatmap.txt")) %>%
  groupOTU(grp)

# Fetch countries
reduced_tree3$tip.country <- get_countries(reduced_tree3, metadata)

# Transform in data.frame
dd2 <- do.call(rbind, Map(data.frame, node=1:length(reduced_tree3$tip.label), Country=reduced_tree3$tip.country))

# Sort countries
dd2$Country <- as.character(dd2$Country)
dd2 <- dd2[order(dd2$Country),]
dd2$Country <- as.factor(dd2$Country)

# Remove row.names
row.names(dd2) <- NULL

# Check results
custom_palette <- 
  createPalette(N = 33,
                seedcolors = c('#b22222', '#7fc97f','#beaed4',
                               '#fdc086','#ffff99','#386cb0',
                               '#f0027f','#bf5b17','#666666'))
p <- reduced_tree3 %>%
  ggtree(layout='rectangular', branch.length = 'none') %<+% dd2 +
  geom_tippoint(aes(color=Country), size=1.2, alpha=1.5) +
  scale_color_manual(values = c(unname(custom_palette), 'black')) +
  geom_tiplab(
    size = 4,
    linesize = .5,
    offset = 0.2,
    align = TRUE,
    aes(subset=label %in% c("Kv15", "Kv35", "Kv57", "Kv97", "Kv104")))
p

# Change label from number to category
gene_matrix_t[gene_matrix_t==0] <- "absent"
gene_matrix_t[gene_matrix_t==1] <- "present"

# Change colors by gene class
classes <- gene_matrix %>% select(Gene_Class) %>% unique()
for (i in as.list(classes$Gene_Class)) {
  genes <- gene_matrix %>% filter(Gene_Class == i) %>% select(Gene) %>% unique()
  
  for (u in as.list(genes$Gene)) {
    gene_matrix_t[gene_matrix_t[, u] == "present", u] <- 
      gsub("present", i, gene_matrix_t[gene_matrix_t[, u] == "present", u])
  }
  
}

# Plot with heatmap
heat_p <- gheatmap(p, gene_matrix_t, colnames_angle = 45, font.size = 4.5, offset = 10, colnames_offset_y = 0.5, 
                   colnames_position = "top", color = "black", hjust = 0, width = 3.5) +
  scale_fill_manual(values=c('white', 'lightgray', '#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69'),
                    name   = "Virulence genes") +
  theme(
    legend.background=element_rect(fill=NA),
    legend.title=element_text(size=13), # The size should be adjusted with different devout.
    legend.text=element_text(size=11),
    legend.spacing.y = unit(0.25, "cm")) +
  ylim(0,185) +
  guides(colour = guide_legend(override.aes = list(size=3))) # adjust size of dots in legend
heat_p

# Save plots
ggsave(filename = "./ggtree_virulence_genes_colored_with_countries.png", plot = heat_p, device = "png")
ggsave(filename = "./ggtree_virulence_genes_colored_with_countries.svg", plot = heat_p, device = "svg", 
       height = 15, width = 20)
