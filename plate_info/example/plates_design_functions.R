# All parsing analysis and plotting functions

# Libraries
library(tidyverse)

# Plates planning functions ----



# Function to distribute strains on new 96-well plates and to take account of strain groups so that strains from the same group are on the same plate
plate_distributor <- function(df_isolates,
                              n_replicates_controls,
                              exclude_outer_wells = F,
                              plate_number_format = NULL,
                              plate_number_start = 1, # numeric with number of first plate (can be used to add new plates to existing collection)
                              plot = T){
  
  
  # check variables
  df <- df_isolates
  if (!"strain_group" %in% colnames(df)){
    df <- df %>%
      mutate(strain_group = NA_character_)
  }
  
  n_isolates_plate <- (96 - sum(n_replicates_controls))/3
  
  df <- df %>%
    arrange(strain_group, sample_id) 
  
  plates_vector <- c() # empty vector with plate numbers. 
  # The for loop will generate three variables
  # x (plate number), j (number of strains on the plate including the current one), 
  # k (number of strains on the plate if all samples from the same patient are added)
  x <- 1
  j <- 1
  for (i in 1:nrow(df)){
    if (i > 1){ # need to start from the second row to avoid NA
      episode <- df$strain_group[i]
      previous <- df$strain_group[i -1]
      
      if (episode != previous ){ # compute k only for the first strain if multiple strains per episode
        n_isolates_episode <- length(which(df$strain_group == df$strain_group[i]))
        k <- j + n_isolates_episode - 1 # substract 1 because first isolates already included in j
        # if the number of strains (including expected multiple same-patient isolates) exceeds the maximum number, start a new plate
        if (j > n_isolates_plate | k > n_isolates_plate){
          j <- 1
          x <- x + 1
        }
      } else {
        if (j > n_isolates_plate){
          j <- 1
          x <- x + 1
        }
      }
    }
    plates_vector[i] <- x
    j <- j + 1
  }
  df_strains_plate <- df %>%
    mutate(plate = plates_vector + plate_number_start - 1) 
  
  if (!is.null(plate_number_format)){
    df_strains_plate <- df_strains_plate %>%
      mutate(plate = str_c(plate_number_format, plate))
  }
    
  
  # Quick check that the same-pair strains are on the same plate
  df_strains_plate %>%
    ggplot(aes(x = fct_rev(strain_group), 
               y = fct_relevel(plate, c("GP8", "GP9")))) +
    geom_tile(fill = "navy") +
    coord_flip() +
    labs(x = "", y = "") +
    theme_bw()
  
  return(df_strains_plate)
}

plate_planner <- function(isolates, # character vector of samples/isolates names
                          n_replicates = 3,
                          controls = c("JE2", "TOX-4", "Complete lysis", "Non infected"), # character vector of control names
                          n_replicates_controls = c(3, 3, 3, 3), # numeric with number of replicates for the controls. Must have same length and order like the controls vector 
                          exclude_outer_wells = FALSE, # should the outer wells be left blank (in case of evaporation issues)?
                          randomise = TRUE, # should the position on the plate be randomised?
                          fixed = NULL) {# character vector of isolates / controls with a fixed position
  
  # Add controls to the isolates vector
  samples <- c(isolates, controls)
  # A vector with the replicates
  replicates <- 1:n_replicates
  # Write vector with samples, repeated three times
  sample_id <- sort(rep(samples, length(replicates)))
  # Set vector length to 96 or 60, depending on whether the outer wells are excluded or not
  # (if there are less samples/ replicates, it will generate NA values,
  # if there are more than  n_samples, the exceeding samples will be excluded)
  if (exclude_outer_wells) {
    max_length <- 60
  } else {
    max_length <- 96
  }
  length(sample_id) <- max_length
  # Write vector with replicate number, repeated for each sample, set length to 96
  replicate <- rep(replicates, length(sample_id))
  length(replicate) <- max_length
  
  # Generate dataframe
  replicates_df <- dplyr::tibble(sample_id, replicate) %>%
    replace_na(list(sample_id = "Non infected", replicate = "Non infected"))
  
  # Optional randomisation
  if (!is.null(fixed)) {
    fixed_df <- replicates_df %>%
      filter(sample_id == fixed)
    replicates_df <- replicates_df %>%
      filter(sample_id != fixed)
  } else {
    fixed_df <- NULL
  }
  if (randomise) {
    set.seed(1)
    replicates_df <- sample_frac(replicates_df)
    replicates_df <- replicates_df %>%
      bind_rows(fixed_df)
  }
  print(replicates_df)
  
  ### Add well in format "A01, A02, ..."
  well <- paste0(sort(rep(LETTERS[1:8], 12)), rep(1:12, 8))
  well <- gsub("^([A-Z])([0-9]){1}$", "\\10\\2", well)
  
  
  ### Unused wells
  if (exclude_outer_wells) {
    well <- str_subset(well,"01|12|A|H", negate = TRUE) 
  } 
  
  # Final dataframe in long format
  well_info <- tibble(well) %>%
    dplyr::bind_cols(replicates_df)
  
  # Fix number of replicates for controls
  # Need to write a better code for this!!
  length(n_replicates_controls) <- 4
  if (n_replicates_controls == c(4, 2, 3, 3)){
    well_info <- well_info %>%
      mutate(replicate = if_else(sample_id == "TOX-4" & replicate == 3,
                                 "4", replicate),
             sample_id = if_else(sample_id == "TOX-4" & replicate == 4, 
                                 "JE2", sample_id),
             )
  }
  
  return(well_info)
}

plate_mapper <- function(well_info){
  plate_info <- well_info %>%
    tidyr::replace_na( replace = list( sample_id = "" ) ) %>%
    tidyr::separate( col = well, into = c( "row", "column" ), sep = 1 ) %>%
    dplyr::select( row, column, sample_id ) %>%
    tidyr::spread( key = column, value = sample_id )
  return(plate_info)
}
