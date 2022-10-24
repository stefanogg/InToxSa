library(tidyverse)
#library(umap)
#library(M3C)
library(uwot)

setwd("Operetta/")

Cell_count_clean.df <- readRDS("Cell_count_clean.df")

unique(Cell_count_clean.df$sample_id)

# select samples to represent in UMAP
Cell_count_UMAP_all.df <- Cell_count_clean.df %>%
  mutate(cell_id = str_c(sample_id, plate_number, plate_replicate, replicate, timepoint, well, Field, Object_No, sep = "_")) %>%
  #filter(Number_of_Spots >0) %>%
  filter(sample_id %in% c("NE581","NE119","NE1532","JE2" )) 

# generate umap data frame
Cell_count_UMAP_data.df <- Cell_count_UMAP_all.df %>%
  #head(1000) %>%
  select(cell_id, Total_Spot_Area, 
         Relative_Spot_Intensity, 
         Number_of_Spots, 
         #Number_of_Spots_per_Area_of_Cell, 
         #Spots_per_Cell_Mean, 
         #`Spots_per_Cell_Morphology_Area_[µm²]`, 
         Spots_per_Cell_Morphology_Roundness) %>%
  tibble::column_to_rownames(var = "cell_id")

# generate umap metadata dataframe
UMAP_metadata.df <- Cell_count_UMAP_all.df %>%
  select(cell_id,
         timepoint,
         plate_replicate,
         sample_id,
         replicate,
         strain_group)

# umap
umap.res <-umap(Cell_count_UMAP_data.df, scale = T,init = "laplacian", spread = .5,
                  
                y = UMAP_metadata.df$sample_id %>% as.factor())

umap.df <- umap.res %>% as.data.frame()
  
ggplot(umap.df, aes(x=V1, y=V2, colour = UMAP_metadata.df$sample_id)) +
  geom_point(size = .3)

ggplot(Cell_count_UMAP_all.df %>% 
         filter(plate_replicate == "n1") %>%
         filter(replicate == 1)
       , aes(x = X, y = Y, colour = Number_of_Spots  , size = Number_of_Spots)) +
  geom_point() +
  scale_size_continuous(range = c(0.1, 5)) +
  scale_colour_distiller(palette = "Spectral") +
  facet_wrap(sample_id ~ timepoint, ncol = 2)
  
