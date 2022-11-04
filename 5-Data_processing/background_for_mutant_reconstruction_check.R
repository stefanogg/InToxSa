library(tidyverse)

PI <- readRDS("Genetic_pairs_analysis_Oct_2020/processed_data/PI/dataframes/parameters/PI_sample_parameters.Rda")

pi <- readRDS("Genetic_pairs_analysis_Oct_2020/processed_data/PI/dataframes/parameters/PI_parameters.Rda")

strain_meta <- read.csv("Genetic_pairs_analysis_Oct_2020/metadata/strain_metadata_corrected_with_clade_CC.csv")

PI_meta <- merge(PI, strain_meta, by = "sample_id")

hist(PI_meta$AUC_death_IQR, breaks = 100)

PI_meta %>% filter(AUC_death_IQR < 50) %>% 
  ggplot2::ggplot(aes(x=as.factor(CC), y= AUC_death_mean, color = as.factor(CC))) +
  geom_boxplot() +
  geom_jitter()



select <- PI_meta %>%
  filter( CC == "239") %>%
select(sample_id, mean_AUC_death, CC, ST, AUC_death_IQR)

# Use BPH3370
# Growth curve and PI curve ok (high PI (~JE2 level) average growth AUC, low variance for botn parameters)