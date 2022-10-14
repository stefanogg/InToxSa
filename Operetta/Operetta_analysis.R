library(tidyverse)

setwd("/home/rguerillot/Documents/Travail/Abdou_project/Staph_infection_project/github_analysis/VANANZ_phenotypes")

# import cell count operetta data
read_tsv_filename <- function(flnm) {
  read_tsv(flnm, skip = 9, trim_ws = T) %>% 
    mutate(filename = flnm)
}

Cell_count.df <-list.files(path = "Operetta/raw_data", pattern = "*Cells.txt", full.names = T, recursive = T) %>% 
  map_df(~read_tsv_filename(.)) 


# Tidy Operetta combined dataframe

names(Cell_count.df) <- gsub(" ", "_", names(Cell_count.df))
names(Cell_count.df) <- gsub("Cells_-_", "", names(Cell_count.df))
Cell_count_clean.df <- Cell_count.df %>%
  separate(col = filename, sep = "/", into = c("a","b","experiment_id","d")) %>%
  separate(col = experiment_id, sep = " ", into = c("plate_number", "cell_type", "timepoint", "plate_replicate"), remove = F) %>%
  select(-Vananz_GP2_3H_n1, -HeLa, -`15000_cells`, -`3h_infection`, -a, -b, -d, -Compound, -Concentration, -Cell_Count, -Cell_Type, -X21) %>%
  mutate(Row = as.character(Row))

raw_to_ABC.df <- data_frame(Row = c(1,2,3,4,5,6,7,8), row = c("A","B","C","D","E","F","G","H")) %>%
  mutate(Row = as.character(Row))

Cell_count_clean.df <- left_join(Cell_count_clean.df, raw_to_ABC.df) %>%
  select(-Row) %>%
  select(Row = row, Column, everything()) %>%
  mutate(Column = str_pad(Column, 2, pad = "0")) %>%
  mutate(well = paste0(Row, Column)) %>%
  select(well, everything())
  
# import plate metadata and merge
plate_info_GP.df <- read.csv("Ideas_Grant_2020_analysis/Genetic_pairs_analysis_completed/plate_info/well_info.csv") %>%
  rename(strain_group = plan) %>%
  mutate(strain_group = ifelse(grepl("BPH", sample_id), "CLINICAL", "CONTROL"))
plate_info_GPV1.df <- read.csv("plate_info/Nebraska_2021/well_info/well_info_plateGPV1.csv")
plate_info.df <- rbind(plate_info_GP.df, plate_info_GPV1.df) %>%
  mutate(replicate = as.factor(replicate))

Cell_count_clean.df <- merge(Cell_count_clean.df , plate_info.df, by = c("plate_number", "well"), all.x = T)

Cell_count_clean.df <- Cell_count_clean.df %>%
  mutate(sample_id = ifelse(sample_id == "Complete lysis", yes = "Non infected", no = as.character(sample_id)))


#saveRDS(Cell_count_clean.df, file = "Operetta/Cell_count_clean.df")


# count cells per strain
t <- Cell_count_clean.df  %>% group_by(experiment_id, sample_id) %>%
  count()

# estimate size of bacteria spot
hist(Cell_count_clean.df$Total_Spot_Area/Cell_count_clean.df$Number_of_Spots, breaks = 100)
hist(Cell_count_clean.df$Total_Spot_Area/Cell_count_clean.df$Number_of_Spots, breaks = 1000, xlim = c(0,40))
# bacteria/spot are = 22

# Count number of bacteria based or spot area
Cell_count_clean.df <- Cell_count_clean.df %>%
  mutate(nb_bacteria = round(Total_Spot_Area/22))

# plot all cell count and bacteria per cells for every experiment ----
library(ggsci)
library(ggridges)

# o cell count infected/non-infected + violin estimated bacteria/cell (spot area based) ----
for (exp in unique(Cell_count_clean.df$experiment_id)) {
  
  df <- Cell_count_clean.df %>%
    filter(experiment_id == exp)  
  
  p1 <-  ggplot(df %>% filter(nb_bacteria >=1), aes(y = sample_id, x= nb_bacteria, fill = "red")) +
    geom_violin(scale = "count", bw = .3) +
    ylab("") +
    xlab("bacteria per infected cells") +
    theme_bw()+
    xlim(1,25)+
    theme(legend.position = "none")
  p1
  
  p2 <- ggplot(df , aes(x = sample_id, fill = nb_bacteria >=1 )) +
    geom_bar(stat = "count", position = "stack") +
    xlab("") +
    ylab("number of cells") +
    ggtitle(exp) +
    theme(plot.title = element_text(size = 15)) +
    scale_fill_manual(name = "Infected:", labels = c("no" , "yes"), values = c("#4575b4", "#d73027")) +
    coord_flip() +
    theme_bw()
  p2
  
  
  p21 <- cowplot::plot_grid(p2,p1, align = "vh", ncol = 2)
  print(p21)
  ggsave(path = "Operetta/plots", filename = paste0(exp, "_violin_estimated"), plot = p21, device = "jpeg", scale = 2.5)
}

# o 3h/24h comparison  cell count infected/non-infected + violin bacteria/cell ----

Cell_count_clean.df <- Cell_count_clean.df %>%
  mutate(plate_replicate_id = paste0(plate_number, "_", plate_replicate))

for (exp in unique(Cell_count_clean.df$plate_replicate_id)) {
  
  df <- Cell_count_clean.df %>%
    filter(plate_replicate_id == exp)
  
  for (sample in unique(df %>% filter(strain_group != "CONTROL") %>% .$sample_id)) {
    
    df1 <- df %>%
      filter(sample_id == sample | strain_group == "CONTROL")
    
    p1 <-  ggplot(df1 %>% filter(nb_bacteria >= 1), aes(y = sample_id, x= nb_bacteria, fill = timepoint)) +
      geom_violin(scale = "count", bw = .3) +
      ylab("") +
      xlab("bacteria per infected cells") +
      theme_bw()+
      theme(legend.position = "none") +
      scale_fill_manual(
        #labels = c("3h", "24h"),
        values = c("#4daf4a", "#ff7f00"),
        name = "Time post infection", guide = "legend") +
      facet_wrap(~sample_id, ncol = 1,strip.position = "left", scales = "free_y") +
      theme_bw()+
      theme(plot.title = element_text(size = 15),
            strip.text.y.left = element_blank(),
            strip.background = element_blank(),
            strip.placement = "outside",
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position = "none"
      ) +
      xlim(1,30)
    p1
    
    p2 <- ggplot(df1 , aes(x = timepoint, fill = nb_bacteria >=1 )) +
      geom_bar(stat = "count", position = "stack") +
      xlab("") +
      ylab("number of cells") +
      ggtitle(exp) +
      theme(plot.title = element_text(size = 15),
            aspect.ratio = 1,
            strip.background = element_blank(),
            strip.placement = "outside"
      ) +
      scale_fill_manual(name = "Infected:", labels = c("no" , "yes"), values = c("#4575b4", "#d73027")) +
      coord_flip() +
      #scale_y_continuous(trans='log2')+
      facet_wrap(~sample_id, ncol = 1,strip.position = "left") +
      theme_bw()+
      theme(plot.title = element_text(size = 15),
            strip.text.y.left = element_text(angle = 0),
            strip.background = element_blank(),
            strip.placement = "outside",
            legend.position = "none") 
    
    
    
    
    p21 <- cowplot::plot_grid(p2,p1, align = "vh", ncol = 2)
    # print(p21)
    ggsave(path = "Operetta/plots", filename = paste0(sample, "_", exp, "_3h_24h_violins_estimated_nb_bacteria_replicates_combined"), plot = p21, device = "jpeg", width = 10, height = 5, scale = 1)
  }
}


# o cell count infected/non-infected + violin bacteria/cell ----
for (exp in unique(Cell_count_clean.df$experiment_id)) {
  
  df <- Cell_count_clean.df %>%
    filter(experiment_id == exp)  
  
p1 <-  ggplot(df, aes(y = sample_id, x= Number_of_Spots, fill = "red")) +
  geom_violin(scale = "width", adjust = .5) +
  ylab("") +
  xlab("bacteria per infected cells") +
  theme_bw()+
  theme(legend.position = "none")
p1

p2 <- ggplot(df , aes(x = sample_id, fill = Number_of_Spots >0 )) +
  geom_bar(stat = "count", position = "stack") +
  xlab("") +
  ylab("number of cells") +
  ggtitle(exp) +
  theme(plot.title = element_text(size = 15)) +
  scale_fill_manual(name = "Infected:", labels = c("no" , "yes"), values = c("#4575b4", "#d73027")) +
  coord_flip() +
  theme_bw()
p2


p21 <- cowplot::plot_grid(p2,p1, align = "vh", ncol = 2)
#print(p21)
ggsave(path = "Operetta/plots", filename = paste0(exp, "_violin_replicates_combined"), plot = p21, device = "jpeg", scale = 2.5)
}

# o 3h/24h comparison  cell count infected/non-infected + violin bacteria/cell ----

Cell_count_clean.df <- Cell_count_clean.df %>%
  mutate(plate_replicate_id = paste0(plate_number, "_", plate_replicate))

for (exp in unique(Cell_count_clean.df$plate_replicate_id)) {
 
  df <- Cell_count_clean.df %>%
    filter(plate_replicate_id == exp)
 
   for (sample in unique(df %>% filter(strain_group != "CONTROL") %>% .$sample_id)) {
     
     df1 <- df %>%
        filter(sample_id == sample | strain_group == "CONTROL")
     
     p1 <-  ggplot(df1 %>% filter(Number_of_Spots > 0), aes(y = sample_id, x= Number_of_Spots, fill = timepoint)) +
       geom_violin(scale = "area", bw = .3) +
       ylab("") +
       xlab("bacteria per infected cells") +
       theme_bw()+
       theme(legend.position = "none") +
       scale_fill_manual(
         #labels = c("3h", "24h"),
         values = c("#4daf4a", "#ff7f00"),
         name = "Time post infection", guide = "legend") +
       facet_wrap(~sample_id, ncol = 1,strip.position = "left", scales = "free_y") +
       theme_bw()+
       theme(plot.title = element_text(size = 15),
             strip.text.y.left = element_blank(),
             strip.background = element_blank(),
             strip.placement = "outside",
             axis.title.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),
             legend.position = "none"
       ) +
       xlim(1,30)
     p1
  
     p2 <- ggplot(df1 , aes(x = timepoint, fill = Number_of_Spots >0 )) +
       geom_bar(stat = "count", position = "stack") +
       xlab("") +
       ylab("number of cells") +
       ggtitle(exp) +
       theme(plot.title = element_text(size = 15),
             aspect.ratio = 1,
             strip.background = element_blank(),
             strip.placement = "outside"
       ) +
       scale_fill_manual(name = "Infected:", labels = c("no" , "yes"), values = c("#4575b4", "#d73027")) +
       coord_flip() +
       #scale_y_continuous(trans='log2')+
       facet_wrap(~sample_id, ncol = 1,strip.position = "left") +
       theme_bw()+
       theme(plot.title = element_text(size = 15),
             strip.text.y.left = element_text(angle = 0),
             strip.background = element_blank(),
             strip.placement = "outside",
             legend.position = "none") 
     
     
     
     
     p21 <- cowplot::plot_grid(p2,p1, align = "vh", ncol = 2)
     # print(p21)
     ggsave(path = "Operetta/plots", filename = paste0(sample, "_", exp, "_3h_24h_violins_replicates_combined"), plot = p21, device = "jpeg", width = 10, height = 5, scale = 1)
   }
}




# o cell count infected/non-infected + ridge bacteria/cell ----
for (exp in unique(Cell_count_clean.df$experiment_id)) {
  
  df <- Cell_count_clean.df %>%
    filter(experiment_id == exp)
  
  p1 <-  ggplot(df) +
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, aes(fill = 0.5- abs(0.5 - stat(ecdf)),  y = sample_id, x = Number_of_Spots), 
                        bandwidth = .5, scale = 4, size = 0.25, alpha = .2, color = "white") +
    ylab("") +
    xlab("bacteria per infected cells") +
    theme_bw()+
    scale_fill_gradient( low = "yellow", high = "#d73027",  name = "probability") +
    #scale_fill_viridis_c(name = "probability", direction = -1) +
    xlim(0, 25)
  p1
  
 
  
  
  p2 <- ggplot(df , aes(x = sample_id, fill = Number_of_Spots >0 )) +
    geom_bar(stat = "count", position = "stack") +
    xlab("") +
    ylab("number of cells") +
    ggtitle(exp) +
    theme(plot.title = element_text(size = 15)) +
    scale_fill_manual(name = "Infected:", labels = c("no" , "yes"), values = c("#4575b4", "#d73027")) +
    coord_flip() +
    theme_bw()
  p2
  
  
  p21 <- cowplot::plot_grid(p2,p1, align = "vh", ncol = 2)
  #print(p21)
  ggsave(path = "Operetta/plots", filename = paste0(exp, "_ridges_replicates_combined"), plot = p21, device = "jpeg", scale = 2.5)
}

# o 3h/24h comparison  cell count infected/non-infected + ridge bacteria/cell ----

Cell_count_clean.df <- Cell_count_clean.df %>%
  mutate(plate_replicate_id = paste0(plate_number, "_", plate_replicate))

for (exp in unique(Cell_count_clean.df$plate_replicate_id)) {
  
  df <- Cell_count_clean.df %>%
    filter(plate_replicate_id == exp)
  
  p1 <-  ggplot(df) +
    geom_density_ridges(aes(fill = timepoint,  
                            y = sample_id, x = Number_of_Spots), 
                        bandwidth = 0.5,
                        scale = 5, size = 0.25, rel_min_height = 0.00005, alpha = .5, color = "white") +
    ylab("") +
    xlab("bacteria per infected cells") +
    scale_fill_manual(
      labels = c("3h", "24h"),
      values = c("#4daf4a", "#ff7f00"),
      name = "Time post infection", guide = "legend"
    ) +
    facet_wrap(~sample_id, ncol = 1,strip.position = "left", scales = "free_y") +
    theme_bw()+
    theme(plot.title = element_text(size = 15),
          strip.text.y.left = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
    )+ 
    #scale_fill_gradient( low = "yellow", high = "#d73027",  name = "probability") +
    #scale_fill_viridis_c(name = "probability", direction = -1) +
    xlim(0, 25)
  p1

  
  
  p2 <- ggplot(df , aes(x = timepoint, fill = Number_of_Spots >0 )) +
    geom_bar(stat = "count", position = "stack") +
    xlab("") +
    ylab("number of cells") +
    ggtitle(exp) +
    theme(plot.title = element_text(size = 15),
          aspect.ratio = 1,
          strip.background = element_blank(),
          strip.placement = "outside"
          ) +
    scale_fill_manual(name = "Infected:", labels = c("no" , "yes"), values = c("#4575b4", "#d73027")) +
    coord_flip() +
    facet_wrap(~sample_id, ncol = 1,strip.position = "left") +
    theme_bw()+
    theme(plot.title = element_text(size = 15),
          strip.text.y.left = element_text(angle = 0),
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.position = "none"
    ) 
  
  
  
  p21 <- cowplot::plot_grid(p2,p1, align = "vh", ncol = 2)
  print(p21)
  ggsave(path = "Operetta/plots", filename = paste0(exp, "_3h_24h_ridges_replicates_combined"), plot = p21, device = "jpeg", width = 10, height = 20, scale = 1)
}





for (exp in unique(Cell_count_clean.df$experiment_id)) {
  
df <- Cell_count_clean.df %>%
  filter(experiment_id == exp)

p1 <- ggplot(df %>% filter(Number_of_Spots>0), aes(x = sample_id, y = Number_of_Spots, colour = replicate)) +
  #geom_boxplot(outlier.size = 0) +
  geom_pointrange(stat = "summary", position = position_dodge(width = .8), size = .2,
                  fun.min = min,
                  fun.max = max,
                  fun = median) +
  coord_flip() +
  xlab("") +
  ylab("bacteria per infected cells") +
  theme_bw()+
  scale_colour_npg() 
  #geom_point(position = position_jitterdodge(), alpha = .3)
p1

p1b <- ggplot(df %>% filter(Number_of_Spots>0), aes(x = sample_id, y = Number_of_Spots)) +
  #geom_boxplot(outlier.size = 0) +
  geom_pointrange(stat = "summary", position = position_dodge(width = .8), size = .2,
                  fun.min = median,
                  fun.max = median,
                  fun = median) +
  coord_flip() +
  xlab("") +
  ylab("bacteria per infected cells") +
  theme_bw()+
  scale_fill_npg() 
#geom_point(position = position_jitterdodge(), alpha = .3)
p1b

p1c <- ggplot(df %>% filter(Number_of_Spots>0), aes(x = sample_id, y = Number_of_Spots/Total_Spot_Area)) +
  #geom_boxplot(outlier.size = 0) +
  geom_pointrange(stat = "summary", position = position_dodge(width = .8), size = .2,
                  fun.min = min,
                  fun.max = max,
                  fun = median) +
  coord_flip() +
  xlab("") +
  ylab("spot per area of spot (infected cells)") +
  theme_bw()+
  scale_fill_npg() 
#geom_point(position = position_jitterdodge(), alpha = .3)
p1c

library(ggridges)
p1d <- ggplot(df) +
  geom_density_ridges(aes(y = sample_id, x = Number_of_Spots, fill = sample_id),scale = 10, size = 0.25, rel_min_height = 0.00005, alpha = .6, color = "white", from = 0, to = 100) +
  xlab("") +
  ylab("bacteria per infected cells") +
  theme_bw()+
  scale_fill_cyclical(
#    breaks = c("1980 Indy", "1980 Unionist"),
#    labels = c(`1980 Indy` = "Indy", `1980 Unionist` = "Unionist"),
    values = c("#00AFBB", "#E7B800")
#    name = "Option", guide = "legend"
  )
p1d  
  

p1e <-  ggplot(df) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, aes(fill = 0.5- abs(0.5 - stat(ecdf)),  y = sample_id, x = Number_of_Spots), 
                      scale = 5, size = 0.25, rel_min_height = 0.00005, alpha = .2, color = "white") +
  ylab("") +
  xlab("bacteria per infected cells") +
  theme_bw()+
  #scale_fill_gradient2( low = "#053061", mid = "#4daf4a", high = "#67001f", midpoint = 0.25, name = "probability") +
  scale_fill_viridis_c(name = "probability", direction = -1) +
  xlim(0, 25)
p1e


p2 <- ggplot(df %>% group_by(replicate, sample_id) %>%
         count(), aes(x = sample_id, y = n, fill = replicate)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("") +
  ylab("number of cells") +
  ggtitle(exp) +
  theme(legend.position = "none",
        plot.title = element_text(size = 15)) +
  scale_fill_npg() +
  coord_flip() 
p2

p3 <- ggplot(df %>% 
               filter(Number_of_Spots >0) %>%
               group_by(replicate, sample_id) %>%
               count(), aes(x = sample_id, y = n, fill = replicate)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("") +
  ylab("number of infected cells") +
  theme(legend.position = "none")+
  scale_fill_npg() +
  coord_flip() 
p3

p123 <- cowplot::plot_grid(p2,p3,p1e, align = "h", rel_widths = c(0.45, 0.45, 0.55), ncol = 3)
print(p123)
ggsave(path = "Operetta/plots", filename = paste0(exp, "_spot_per_area_of_spot"), plot = p123, device = "jpeg", scale = 2.5)
}
