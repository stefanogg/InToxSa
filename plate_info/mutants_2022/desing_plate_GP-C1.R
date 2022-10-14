# Here we design plate GP-C1

# We start with a dataframe of isolates (that may be grouped), a list of controls. 
# The dataframe is a adapted from Abdou's email

# Library
library(tidyverse)

# free environment
rm(list = ls())

setwd("plate_info/mutants_2022/")

# Source functions
source("plates_design_functions.R")

# Set working directory
setwd(str_c(here::here(), "/plate_info/mutants_2022"))

# Raw data and variables
df_isolates <- read_csv("df_isolates.csv")
controls <- c("JE2", "TOX-4", "BPH3370", "Complete lysis", "Non infected")
n_replicates_controls <- c(5, 5, 5, 3, 3) # a numeric of the same length as the control vector for the number of replicates for each control

# Distribute strains across plates. This function returns a dataframe of the same length of df_isolates
# with a supplementary column indicating the plate assignment (as plate number). If plot = T (the default) it prints a plot to check the distribution
df_plate_distribution <- plate_distributor(df_isolates = df_isolates,
                  n_replicates_controls = n_replicates_controls,
                  plate_number_format = "GP-C")

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
    n_replicates_controls = n_replicates_controls,
    exclude_outer_wells = T,
  ) %>%
    mutate(plate_number = x)
  
  write_csv(well_info,
            str_c(dir, "well_info_plate", x, ".csv"))
  
  return(well_info)
}) %>%
  bind_rows()

# Fix replicate number for controls
set.seed(1234)
well_info <- well_info %>%
  slice_sample(prop = 1) %>%
  group_by(sample_id) %>%
  mutate(replicate = if_else(sample_id %in% controls, row_number(), replicate)) %>%
  arrange(well)

well_info %>%
  write_csv("well_info/well_info_plateGP-C1.csv")

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

