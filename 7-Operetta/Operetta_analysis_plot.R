library(tidyverse)

# set working directory
setwd(paste0(here(),"/7-Operetta/"))

# import highcontent imaging data
Cell_count_corrected.df <- readRDS(file = "Cell_count_corrected.df.Rds") %>%
  mutate(cell_id = str_c(sample_id, plate_number, plate_replicate, replicate, timepoint, well, Field, Object_No, sep = "_")) %>%
  mutate(timepoint = factor(.$timepoint, levels = c("3h", "24h"))) %>%
  mutate(plate_replicate_id = paste0(plate_number, "_", plate_replicate)) 

# plot %cell infected and bacteria/cell heatmap

df2 <- Cell_count_corrected.df %>%
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

df3 <- Cell_count_corrected.df %>%
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


p23 <- cowplot::plot_grid(p2,p3, align = "vh", ncol = 2)
p23
#ggsave(path = "plots", filename = paste0(exp, "_heatmap.eps"), plot = p21, device = "eps" , width = 6.5, height = 6)

