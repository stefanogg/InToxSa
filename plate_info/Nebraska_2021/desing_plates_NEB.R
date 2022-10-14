# Design NEB plates based on convergence analysis

# We start with a dataframe of isolates (that may be grouped), a list of controls. We assume 3 replicates for the strains, but the number of replicates 
# of controls can vary

# Library
library(tidyverse)
library(readxl)
rm(list = ls())

# Set working directory
setwd(str_c(here::here(), "/plate_info/Nebraska_2021/"))

# Source functions
source("../../Functions/all_functions.R")

# Raw data and variables
df_isolates <- read_xlsx("../../Genetic_pairs_analysis_Oct_2020/processed_data/convergence/selection_convergent_genes_neb_df.xlsx", sheet = 2) %>%
  transmute(sample_id = neb_mutant_id,
            strain_group = "NEBRASKA")
controls <- c("JE2", "Complete lysis", "Non infected")
n_replicates_controls <- c(3, 3, 3) # a numeric of the same length as the control vector for the number of replicates for each control

# Distribute strains across plates. This function returns a dataframe of the same length of df_isolates
# with a supplementary column indicating the plate assignment (as plate number). If plot = T (the default) it prints a plot to check the distribution
df_plate_distribution <- plate_distributor(df_isolates = df_isolates,
                  n_replicates_controls = n_replicates_controls,
                  plate_number_format = "GPV") 

# Randomise strains on the plate. This functions returns a dataframe of length df_isolates*3 with the well_id of each strain replicate.
# It also saves individual well_info dataframes

dir <- "well_info/"
dir.create(dir)

well_info <- lapply(unique(df_plate_distribution$plate), function(x){
  df <- df_plate_distribution %>%
    filter(plate == x)
  well_info <- plate_planner(
    isolates = df$sample_id,
    randomise = TRUE,
    fixed = "Complete lysis",
    controls = controls,
    n_replicates_controls = n_replicates_controls
  ) %>%
    mutate(plate_number = x)
  
  # add strain group 
  well_info <- well_info %>%
    left_join(df_plate_distribution %>% select(-plate)) %>%
    mutate(strain_group = if_else(sample_id %in% controls,
                                  "CONTROL",
                                  strain_group))
  
  write_csv(well_info,
            str_c(dir, "well_info_plate", x, ".csv"))
  
  return(well_info)
}) %>%
  bind_rows()

# Manually correct two wells (see Abdou on teams, 18/03/2021)
well_info <- well_info %>%
  mutate(sample_id = case_when(
    well == "A09" ~ "NE369",
    well == "E12" ~ "NE362",
    TRUE ~ sample_id
  ))

well_info %>%
  write_csv("well_info/well_info_plateGPV1.csv")

# Generate plate maps in 96-well format
dir <- "plate_maps/"
dir.create(dir)

plate_info <- lapply(unique(df_plate_distribution$plate), function(x){
  df <- well_info %>%
    filter(plate_number == x)
  plate_info <- plate_mapper(df)
  return(plate_info)
})

names(plate_info) <- str_c(
  "plate_", 
  unique(df_plate_distribution$plate))


purrr::walk(names(plate_info), function(x){
  df <- plate_info[[x]]
  file <- str_c(
    dir, "/",
    x,
    ".csv"
  )
  write_csv(df, file)
})

