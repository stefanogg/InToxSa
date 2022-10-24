
library(tidyverse)

setwd("/home/rguerillot/Documents/Travail/Abdou_project/Staph_infection_project/github_analysis/VANANZ_phenotypes")

# import cell count operetta data
read_tsv_filename <- function(flnm) {
  read_tsv(flnm, skip = 9, trim_ws = T) %>% 
    mutate(filename = flnm)
}

Spot_count.df <-list.files(path = "Operetta/raw_data", pattern = "*Spots.txt", full.names = T, recursive = T) %>% 
  map_df(~read_tsv_filename(.)) 


# Tidy Operetta combined dataframe

names(Spot_count.df) <- gsub(" ", "_", names(Spot_count.df))
names(Spot_count.df) <- gsub("Cells_-_", "", names(Spot_count.df))
Spot_count_clean.df <- Spot_count.df %>%
  #head(10000) %>%
  separate(col = filename, sep = "/", into = c("a","b","experiment_id","d")) %>%
  separate(col = experiment_id, sep = " ", into = c("plate_number", "cell_type", "timepoint", "plate_replicate"), remove = F) %>%
  select( -Vananz_GP2_3H_n1, -HeLa, -`15000_cells`, -`3h_infection`, -a, -b, -d, -Compound, -Concentration, -Cell_Count, -Cell_Type, -X23) %>%
  mutate(Row = as.character(Row))

raw_to_ABC.df <- data_frame(Row = c(1,2,3,4,5,6,7,8), row = c("A","B","C","D","E","F","G","H")) %>%
  mutate(Row = as.character(Row))

Spot_count_clean.df <- left_join(Spot_count_clean.df, raw_to_ABC.df) %>%
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

Spot_count_clean.df <- merge(Spot_count_clean.df , plate_info.df, by = c("plate_number", "well"), all.x = T)

Spot_count_clean.df <- Spot_count_clean.df %>%
  mutate(sample_id = ifelse(sample_id == "Complete lysis", yes = "Non infected", no = as.character(sample_id)))



# count cells per strain
t <- Spot_count_clean.df  %>% group_by(experiment_id, sample_id) %>%
  count()

# estimate size of bacteria spot
hist(Spot_count_clean.df$`Spots_-_Spot_Area_[px²]`, breaks = 1000, xlim = c(0,60))

ggplot(Spot_count_clean.df, aes(x = `Spots_-_Spot_Area_[px²]`)) +
  geom_density()+
  facet_wrap(~sample_id)

ggplot(Spot_count_clean.df, aes(x = `Spots_-_Corrected_Spot__Intensity`)) +
  geom_density()+
  facet_wrap(~sample_id)

# Keep spots/bacteria between 15 and 30 only
Spot_count_clean.df2 <- Spot_count_clean.df %>%
  filter(`Spots_-_Spot_Area_[px²]` > 15) %>%
  filter(`Spots_-_Spot_Area_[px²]` < 30) 

t <- count(Spot_count_clean.df2, timepoint, sample_id, plate_replicate, replicate)   

ggplot(Spot_count_clean.df2, aes(x = `Spots_-_Spot_Area_[px²]`)) +
  geom_density() +
  facet_wrap(~sample_id)

ggplot(Spot_count_clean.df2, aes(x = Spot_count_clean.df2$`Spots_-_Uncorrected_Spot__Peak_Intensity`)) +
  geom_density() +
  facet_wrap(~sample_id)

# Create unique cell_id and count spots bacteria for each cells
Spot_Cell_count.df <- Spot_count_clean.df2 %>%
  mutate(cell_id = str_c(sample_id,
                         plate_number, 
                         plate_replicate, 
                         replicate, 
                         timepoint, 
                         well, 
                         Field, 
                         `Spots_-_Object_No_in_Cells`, sep = "_")) %>%
  count(cell_id)

saveRDS(Spot_count_clean.df2, file = "Operetta/Spot_count_clean.df2.Rda")
saveRDS(Spot_Cell_count.df, file = "Operetta/Spot_Cell_count.df.Rda")

