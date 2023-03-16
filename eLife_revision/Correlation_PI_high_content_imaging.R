library(tidyverse)
library(ggsci)
library(ggridges)

setwd(paste0(here(),"/7-Operetta/"))


# Cell_count_corrected.df <- readRDS(file = "Cell_count_corrected.df.Rda") %>%
#   filter(experiment_id %in% c("GPV1 HeLa 3h n1", "GPV1 HeLa 24h n1"))

#saveRDS(Cell_count_corrected.df, "Cell_count_corrected.df.Rds")

Cell_count_corrected.df <- readRDS(file = "Cell_count_corrected.df.Rds") %>%
  mutate(cell_id = str_c(sample_id, plate_number, plate_replicate, replicate, timepoint, well, Field, Object_No, sep = "_")) %>%
  mutate(timepoint = factor(.$timepoint, levels = c("3h", "24h"))) %>%
  mutate(plate_replicate_id = paste0(plate_number, "_", plate_replicate)) 

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

read.csv("neb_tn_insertions.csv



# o cell count infected/non-infected + violin bacteria/cell ----
for (exp in unique(Cell_count_corrected.df$plate_replicate_id)) {
  
  df <- Cell_count_corrected.df %>%
    filter(plate_replicate_id == exp)  %>%
    mutate(sample_id = ifelse(sample_id == "TOX-4", yes = "PAM-agrA", no = sample_id)) %>%
    filter(sample_id != "Non infected")
  
  p1 <-  ggplot(df, aes(y = sample_id, x= Number_of_Spots, fill = "red")) +
    geom_violin(scale = "width", adjust = .5) +
    ylab("") +
    xlab("bacteria per infected cells") +
    theme_bw()+
    theme(legend.position = "none")
  p1
  
  df2 <- df %>%
    mutate(Infected = ifelse(Number_of_Spots>0, yes = "yes", no = "no")) %>%
    group_by(sample_id, timepoint,Infected) %>% 
    summarise(count = n()) %>% 
    mutate(perc = count/sum(count)*100)
  
  p2 <- ggplot(df2 , aes(x = sample_id, y = timepoint,  fill = perc )) +
    geom_tile(stat = "identity") +
    xlab("") +
    ylab("") +
    ggtitle(exp) +
    theme(plot.title = element_text(size = 15)) +
    #scale_fill_manual(name = "Infected:", labels = c("no" , "yes"), values = c("#4575b4", "#d73027")) +
    coord_flip() +
    theme_bw() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.text = element_text(size = 12), 
          plot.title = element_text(size = 16, face = 'bold'), 
          axis.title = element_text(face = 'bold')) + 
    scale_fill_gradient(low = "white", high = "steelblue", name = "% infected")
  p2
  
  df3 <- df %>%
    filter(cor_nb_bacteria > 0) %>%
    group_by(sample_id, timepoint) %>% 
    summarise(bacteria_nb_per_cell = mean(cor_nb_bacteria)) 
  
  p3 <- ggplot(df3 , aes(x = sample_id, y = timepoint,  fill = bacteria_nb_per_cell )) +
    geom_tile(stat = "identity") +
    xlab("") +
    ylab("") +
   # ggtitle(exp) +
    theme(plot.title = element_text(size = 15)) +
    #scale_fill_manual(name = "Infected:", labels = c("no" , "yes"), values = c("#4575b4", "#d73027")) +
    coord_flip() +
    theme_bw() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.text = element_text(size = 12), 
          plot.title = element_text(size = 16, face = 'bold'), 
          axis.title = element_text(face = 'bold')) + 
    scale_fill_gradient(low = "white", high = "red", name = "SA/infected cell")
  p3
  
  
  p21 <- cowplot::plot_grid(p2,p3, align = "vh", ncol = 2)
  print(p21)
  ggsave(path = "plots", filename = paste0(exp, "_heatmap.eps"), plot = p21, device = "eps" , width = 6.5, height = 6)
}



