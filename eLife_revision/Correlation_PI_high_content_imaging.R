library(tidyverse)
library(ggpubr)

# set working directory
setwd(here())

# import nebraska tn mutant high content imaging results 
Cell_count_corrected.df <- readRDS(file = "7-Operetta/Cell_count_corrected.df.Rds") %>%
  mutate(cell_id = str_c(sample_id, plate_number, plate_replicate, replicate, timepoint, well, Field, Object_No, sep = "_")) %>%
  mutate(timepoint = factor(.$timepoint, levels = c("3h", "24h"))) %>%
  mutate(plate_replicate_id = paste0(plate_number, "_", plate_replicate)) 

# get %cell infected 24h
perc_infected_24h.df <- Cell_count_corrected.df %>%
  #filter(sample_id != "Non infected") %>%
    # set % infected
  mutate(Infected = ifelse(Number_of_Spots>0, yes = "yes", no = "no")) %>%
  group_by(sample_id, timepoint,Infected) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)*100) %>%
  ungroup() %>%
  filter(timepoint == "24h" & Infected == "yes")

# merge with PI AUC
neb_PI.df <- read.csv("7-Operetta/PI_nebraska.csv")
signif_PI_str <- c("NE1532", "NE119", "NE188", "NE873", "NE1140", "NE117")
perc_infected_24h_PI.df <- merge(perc_infected_24h.df, neb_PI.df, by = "sample_id") %>%
  mutate(label = ifelse(sample_id %in% signif_PI_str, paste0(sample_id,"-", unique_gene_symbol), ""))
  
# get mean number bacteria / infected cells

bact_per_infected_cell_24h.df <- Cell_count_corrected.df %>%
  filter(cor_nb_bacteria > 0) %>%
  group_by(sample_id, timepoint) %>%
  summarise(bacteria_nb_per_cell = mean(cor_nb_bacteria)) %>%
  filter(timepoint == "24h")
  
# merge with PI AUC
bact_per_infected_cell_24h_PI.df <- merge(bact_per_infected_cell_24h.df, neb_PI.df, by = "sample_id") %>%
  mutate(label = ifelse(sample_id %in% signif_PI_str, paste0(sample_id,"-", unique_gene_symbol), ""))

# check correlations

p1 <- ggscatter(perc_infected_24h_PI.df,
                x = "AUC_death_mean", 
                y = "perc",
                color = "steelblue",
                #fill = "lightgray",
                conf.int = T,
                add = "reg.line")+
  geom_point() +
  stat_cor(label.x = 140, label.y = 55, label.sep = "\n") +
  xlab("PI uptake (AUC)") +
  ylab("% Cells infected at 24h") +
  ggrepel::geom_label_repel(data = perc_infected_24h_PI.df, aes(label = label))
p1

p2 <- ggscatter(bact_per_infected_cell_24h_PI.df,
                x = "AUC_death_mean", 
                y = "bacteria_nb_per_cell",
                color = "red",
                #fill = "lightgray",
                conf.int = T,
                add = "reg.line")+
  geom_point() +
  stat_cor(label.x = 140, label.y = 2.5, label.sep = "\n") +
  xlab("PI uptake (AUC)") +
  ylab("*S. aureus* infected/cell at 24h ") +
  ggrepel::geom_label_repel(data = bact_per_infected_cell_24h_PI.df, aes(label = label))+
  theme(axis.title.y = ggtext::element_markdown())
p2

p <- cowplot::plot_grid(p1, p2, ncol = 1)
p

ggsave(file = "Figure5_figure_supplement5.pdf",
       path = "8-eLife_revision/Figure5_prepare/figures/",
       device = "pdf", 
       plot = p, 
       width = 5, 
       height = 8)
