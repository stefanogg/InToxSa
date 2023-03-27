# ---
# title: "Z-factor calculation of the InToxSa assay"
# author: "Romain Guerillot"
# date: "10/03/2023"
# ---

library(tidyverse)

# Import data
PI_test_parameters <- readRDS("../../5-Data_processing/processed_data/allelic_exchange/PI_parameters/PI_sample_parameters_plateGP-C1.Rda")

# write Z-factor calculation function that use PI parameters dataframe as input
get_Zfactor <- function(param_df, pos_control_samp, neg_control_samp, signal){
  
  param_df_pos <- param_df %>% filter(sample_id == pos_control_samp) 
  sd_pos = param_df_pos %>% 
    select(ends_with("sd")) %>% 
    select(starts_with(signal)) %>% 
    print(.) %>%
    pull()
  mean_pos = param_df_pos %>% 
    select(ends_with("mean")) %>% 
    select(starts_with(signal)) %>% 
    print(.) %>%
    pull()
  
  param_df_neg <- param_df %>% filter(sample_id == neg_control_samp) 
  sd_neg = param_df_neg %>% 
    select(ends_with("sd")) %>% 
    select(starts_with(signal)) %>% 
    print(.) %>%
    pull()
  mean_neg = param_df_neg %>% 
    select(ends_with("mean")) %>% 
    select(starts_with(signal)) %>% 
    print(.) %>%
    pull()
  
  zfactor = 1 - (3*sd_pos + 3*sd_neg)/(abs(mean_pos - mean_neg))
  print(paste0("Z-factor = ", zfactor))
  return(data.frame(controls = paste0(pos_control_samp, "-", neg_control_samp),
                   signal = signal,
                   zfactor = zfactor))
}

# Calculate Z-factors

zscore.df <- data.frame()
for (pos in c("JE2", "Complete lysis")) {
  for (neg in c("NE1532", "Non infected")) {
    for (sig in c("AUC_death", "max_death", "max_rate_death")) {
      d <- get_Zfactor(param_df = PI_test_parameters,
                       pos_control_samp = pos,
                       neg_control_samp = neg,
                       signal = sig)
      zscore.df <- rbind(zscore.df, d)
    }
  }
}

zscore.df
