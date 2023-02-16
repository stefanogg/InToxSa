# Here we implement a phylogeny-based approach to infer clonal complexes

# This approach works well for large and medium-sized clades, however, it generates several singletons.
# Other issue: it depends on the phylogeny and the datasets, making comparison more difficult

# Alternative approaches: use PopPUNK or bespoke mash-based clustering

# Library
library(tidyverse)
library(ggtree)

# Working directory
setwd(str_c(here::here(), "/Phylogeny"))

# Clear existing data
rm(list = ls())

# Import and process tree
dir <- "raw_data/"
f <- "~/Documents/Transfer_with_server/core90.aln.treefile"
file.copy(f, dir, overwrite = F)
tree <- read.tree(str_c(dir, basename(f)))
tree <- tree %>%
  phytools::midpoint.root() %>%
  ape::drop.tip("Reference")
tree$edge.length <- tree$edge.length*3.0496e6

# Import sequence type
st_data <- read_csv("../plate_info/strain_metadata.csv") %>%
  transmute(sample_id, 
            ST = as.factor(ST)) %>%
  filter(sample_id %in% tree$tip.label)

# Test: define CC239 clade
tips <- st_data %>%
  filter(ST == "239") %>%
  .$sample_id
mrca <- MRCA(tree, tips)

ggtree(tree) +
  geom_point2(aes(subset = node == mrca), colour = "red", size = 5)

descendants <- tidytree::offspring(tree, mrca, tiponly = T)
descendants <- tree$tip.label[descendants]

treeio::tree_subset(tree, node = mrca, levels_back = 0) %>%
  ggtree() %<+% st_data  +
  geom_tiplab(aes(label = ST), size = 3, align = T, linesize = .1)

df <- tibble(sample_id = descendants, CC = "239")

# Iterate over first 10 ST
st <- st_data %>%
  filter(ST != "-") %>%
  count(ST, sort = T) %>%
  # slice_head(n = 10) %>%
  .$ST
st <- as.character(st)

df_cc <- tibble(sample_id = character(), CC = character(), mrca = numeric())

st_info <- st_data

for (i in st){
  tips <- st_info %>%
    filter(ST == i) %>%
    .$sample_id
  
  if (length(tips) == 0) next()
  
  mrca <- MRCA(tree, tips)
  
  # ggtree(tree) +
  #   geom_point2(aes(subset = node == mrca), colour = "red", size = 5)
  
  descendants <- tidytree::offspring(tree, mrca, tiponly = T)
  descendants <- tree$tip.label[descendants]
  
  # treeio::tree_subset(tree, node = mrca, levels_back = 0) %>%
  #   ggtree() %<+% st_data  +
  #   geom_tiplab(aes(label = ST), size = 3, align = T, linesize = .1)
  
  df_cc <- df_cc %>%
    add_row(sample_id = descendants, CC = i, mrca = mrca)
  
  st_info <- st_info %>%
    anti_join(df_cc)
  
}

# Merge with ST dataframe
cc_data <- st_data %>%
  full_join(df_cc) %>%
  mutate(CC = as.factor(CC))

# Plot on tree
ggtree(tree, layout = "circular") %<+% cc_data  +
  geom_tiplab(aes(label = CC), size = 3, align = T, linesize = .0001) +
  geom_tiplab(aes(label = ST), size = 3, align = T, linesize = .0001, offset = 1e4)

# Save processed data
cc_data %>%
  saveRDS("processed_data/cc_phylo_inferred.Rda")
tree %>%
  treeio::write.tree("processed_data/clean.core90.aln.tree")
