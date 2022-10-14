# All parsing analysis and plotting functions

# Libraries
library(tidyverse)
library(magrittr)
library(readxl)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(broom)
library(ggpubr)

# Plates planning functions ----

# Wrapper of plates planning function to create a pipeline. Work in progress
plate_designer <- function(df_isolates,
                           controls,
                           n_replicates_controls,
                           exclude_outer_wells,
                           randomise = T,
                           fixed = NULL,
                           plot = T,
                           save_dir){
  
}

# Function to distribute strains on new 96-well plates and to take account of strain groups so that strains from the same group are on the same plate
# This work in progress!!!
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
                          exclude_outer_wells = FALSE, # numeric with number of outer wells rings to be left blank (in case of evaporation issues or other reason to keep plate small)?
                          randomise = TRUE, # should the position on the plate be randomised?
                          fixed = NULL) {# character vector of isolates / controls with a fixed position
  
  print(exclude_outer_wells)
  # Add controls to the isolates vector
  samples <- c(isolates, controls)
  # A vector with the replicates
  replicates <- 1:n_replicates
  # Write vector with samples, repeated three times
  sample_id <- sort(rep(samples, length(replicates)))
  # Set vector length to 96 or 60, depending on whether the outer wells are excluded or not
  # (if there are less samples/ replicates, it will generate NA values,
  # if there are more than  n_samples, the exceeding samples will be excluded)
  if (exclude_outer_wells != FALSE) {
    if (exclude_outer_wells == 1) max_length <- 60
    if (exclude_outer_wells == 2) max_length <- 32
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
  if (exclude_outer_wells != FALSE) {
   if (exclude_outer_wells == 1) well <- str_subset(well,"01|12|A|H", negate = TRUE) 
   if (exclude_outer_wells == 2) well <- str_subset(well,"01|02|11|12|A|B|G|H", negate = TRUE) 
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


# plate_planner <- function(isolates, # character vector of samples/isolates names
#                           n_replicates = 3,
#                           controls = c("JE2", "TOX-4", "Complete lysis", "Non infected"), # chacter vector of control names
#                           n_replicates_controls = c(3, 3, 3, 3), # numeric with number of replicates for the controls. Must have same length and order like the controls vector 
#                           exclude_outer_wells = FALSE, # should the outer wells be left blank (in case of evaporation issues)?
#                           randomise = TRUE, # should the position on the plate be randomised?
#                           fixed = "Complete lysis") # character vector of isolates / controls with a fixed position
# {    
#     if (is.character(title)) {
#       my_title <- title
#     } else if (combine_experiments == "plate") {
#       my_title <- str_c("Combined Plate ", exp,  " - ",  type)
#     } else {
#       my_plates <- str_c(unique(df_exp$plate), collapse = "+")
#       my_title <- str_c(exp, " Plate ", my_plates, " - ",  type)
#     }
#     
#     p <- df_exp %>%
#       ggplot(aes(x = timepoint/60, y = OD, label = well_id)) +
#       #stat_summary(geom = "ribbon", fun.data = mean_cl_normal, fun.args=list(conf.int=0.95), fill = "#045a8d", alpha = 0.2) + # Manually checked and give same results, but 95% CI appear to be too narrow 
#       stat_summary(geom = 'ribbon', fun.data = 'mean_sdl', fun.args = list(mult = 1), fill = "#045a8d", alpha = 0.2) + # Use standard deviation instead
#       stat_summary(geom = "point", fun = "mean", color = "#045a8d", size = 1) +
#       facet_wrap(~ sample_id, nrow = nrow, ncol = ncol) +
#       labs(title = my_title,
#            x = "Time (hour)",
#            y = "Fluorescence intensity") +
#       theme_bw()
#     
#     if (annotate) {
#       require(ggrepel)
#       p <- p + 
#         geom_text_repel(data = subset(df_exp,
#                                       timepoint == (count(df_exp, timepoint) %>% 
#                                                       filter(n == max(n)) %>% 
#                                                       filter(timepoint == max(timepoint)) %>% .$timepoint)),
#                         force = 1,
#                         segment.size  = 0.5,
#                         segment.color = "black",
#                         direction     = "both")
#     }
#     
#     if (comparator != F) {
#       comparator_df <- df_exp %>%
#         filter(sample_id == comparator) %>%
#         select(-sample_id)
#       p <- p + 
#         stat_summary(data = comparator_df, geom = 'ribbon', fun.data = 'mean_sdl', fun.args = list(mult = 1), fill = "#b30000", alpha = 0.2) +
#         stat_summary(data = comparator_df, geom = "point", fun = "mean", color = "#b30000", alpha = 0.2, size = 1)
#     }
#     
#     if (all_replicates != F) {
#       if (is.numeric(all_replicates)) { # if all_replicate is a numeric set fitering value accordingly
#         filt_div = all_replicates
#       }else{
#         filt_div = 15 # plot measure every 15 min by default (all_replicate = T)
#       }
#       p <- p +
#         geom_point(data = df_exp %>% filter(divisible(timepoint, filt_div)), size = .05, color = "#045a8d")
#       if (comparator != F) {
#         p <- p + 
#           geom_point(data = comparator_df %>% filter(divisible(timepoint, filt_div)), size = .05, color = "#b30000")
#       }
#     }
#     
#     # Save the plot
#     
#     if (save_plot != F) {
#       
#       # check directory
#       if (is.character(save_plot)) {
#         my_dir <- save_plot
#       } else {
#         my_dir <- "figures"
#       }
#       if (!dir.exists(my_dir))
#         dir.create(my_dir) # check if directory already exists (to avoid error)
#       
#       # create file name
#       my_title_file <- str_replace_all(my_title, "\\s+|-", "_") %>%
#         str_replace_all("(?<=_)_", "")
#       file <- str_c(
#         my_dir,
#         "/Growth_curves_",
#         my_title_file,
#         ".pdf"
#       )
#       
#       # check if plot size is given as an argument
#       if (is.null(plot_size)) {
#         my_width <- NA
#         my_height <- NA
#       } else {
#         my_width <- plot_size[1]
#         my_height <- plot_size[2]
#       }
#       
#       # save
#       ggsave(file,
#              width = my_width,
#              height = my_height,
#              limitsize = FALSE)
#     }
#     
#     
#     all_p = c(all_p, p)
#     
#     if (print_plot) print(p)
#     
# 
#   if (return_plots == T) {
#     return(p)
#   }
# }


plate_parser <- function(plate_map){
  names(plate_map)[1] <- "row"
  plate_data <- plate_map %>%
    gather(key = column, value = result, num_range(prefix = "", range=1:12)) %>%
    mutate(column = formatC(as.numeric(column), width = 2, format = "d", flag = "0")) %>%
    unite(col = "well", row, column, sep = "")
  return(plate_data)
}

# Raw data parsing functions ----

# Parse PI kinetics

# Parse inocula OD

parse_inocula_OD <- function(file,
                             n_skip = 13,
                             col_names = 
                               c("well_row", "well_col", "sample", "inoculum_OD"),
                             normalised = FALSE,
                             exp,
                             plate){
  
  # fix arguments for parsing
  dir  <- getwd()
  path <- str_c(dir,
                file)
  
  
  
  # parse data
  data <- read_xlsx(file, skip = n_skip, col_names = col_names) %>%
    mutate(well = str_c(well_row, 
                        formatC(well_col, 
                                width = 2, 
                                format = "d", 
                                flag = "0")),
           inoculum_OD_blank_normalised = normalised,
           experiment = exp,
           plate = plate) %>%
    select(-c(well_row, well_col, sample)) 
  
  return(data)
}

annotate_wells <- function(df, # dataframe to annotate
                           well_info = "plate_info/well_info_n843.csv", # dataframe mapping strains to wells
                           strain_group = "plate_info/strain_metadata.csv", # dataframe with strain group
                           annotation = "plate_info/strain_metadata.csv", # metadata
                           type = c("minimal", "wide")) {
  
  type <- match.arg(type) # use match.arg() to extract the default ("minimal") if no value is given
  my_annot <- annotation

  # retrieve well info and modify
  controls <- c("JE2", "TOX-4", "Complete lysis", "Non infected")

  well_info <- read_csv(well_info) 
  
  if (is.null(strain_group)){
    well_info <-  well_info %>%
      transmute(sample_id, well, plate = as.character(plate_number), strain_group) %>%
      mutate(
        sample_id = if_else(sample_id == "TOX-5", "JE2", sample_id),
        sample_id = if_else(sample_id == "Total LDH", "Complete lysis", sample_id)
      ) %>%
      mutate(strain_group = if_else(sample_id %in% controls, "CONTROL", strain_group))
  } else {
    strain_group <- read_csv(strain_group)
    well_info <-  well_info %>%
      left_join(strain_group) %>%
      transmute(sample_id, well, plate = as.character(plate_number), strain_group) %>%
      mutate(
        sample_id = if_else(sample_id == "TOX-5", "JE2", sample_id),
        sample_id = if_else(sample_id == "Total LDH", "Complete lysis", sample_id)
      ) %>%
      mutate(strain_group = if_else(sample_id %in% controls, "CONTROL", strain_group))
    
  }

    
  # retrieve metadata
  if (my_annot == "none"){
    df_out <- df %>% left_join(well_info) 
  } else {
    strain_info <- read_csv(annotation) 
    
    if (type == "minimal") strain_info <- select(strain_info, sample_id, strain_group)
    
    # merge
    df_out <- df %>% left_join(well_info) %>% left_join(strain_info)
  }
  
  
  
  return(df_out)
  
  
}



# quick plotting function to check raw data
quick_plot  <- function(df){
  p <- df %>%
    ggplot(aes(x = timepoint, y = f_signal)) +
    geom_point(size = .5, color = "blue") +
    facet_grid(well_row ~ well_col) +
    theme_bw()
  return(p)
}

# parsing function 
parse_kinetics <- function(file, # path to Excel file
                           n_skip = 13, # number of lines to skip (includes Excel column names)
                           n_col = NULL, # numeric with the number of columns to parse
                           interval = 6, # numeric with interval between measurements (in minutes)
                           duration = NULL, # numeric vector with duration of 
                           # the experiment in hours, minutes (e.g. 22 hours and 12 minutes: 22, 12)
                           n_cycles = NULL, # numeric with number of measurements cycles (alternative to duration)
                           wavelength = c(493, 535), # vector of the wavelengths in raw data (numeric or vector of numeric)  
                           keep_wavelength = c(535), # wavelengths to keep (numeric or vector of numeric)
                           exp,
                           cell_type = "HeLa",
                           cell_number = 4e4,
                           plate,
                           plot_data = TRUE,
                           annotate = c("minimal", "wide", "none")){
  # would be better to change argument name as annotate option means something different in the plotting function
  
  # Starting messages
  print(glue::glue("Parsing PI kinetics data from file {file}"))
  
  # Make annotate an optional argument
  if(!missing(annotate)) {
    annotate <- match.arg(annotate)
  }else{
    annotate <- NULL
  }
  
  # Timepoint calculation in minutes
  if (!is.null(duration)) {
    duration_h <- duration[1]*60 # hour component of duration, transformed in minutes
    duration_min <- duration[2] # minute component of duration
    timepoint_max <- duration_h + duration_min
    # generates vector of timepoints (in minutes) of fluorescence measurement. Here: every 6 minutes until 22 hours and 18 minutes
    timepoints <- seq(from = 0, to = timepoint_max, by = interval) 
  } else if (!is.null(n_cycles)) {
    timepoints <- seq(from = 0, by = interval, length.out = n_cycles) 
  }
  
  
  # generate vector of column names based on timepoints and wavelength. Here, both wave length were used
  col_names_WL <- sapply(wavelength, function(x){
    str_c("WL", x, "_", timepoints)
  })
  col_names <- c(
    "well",
    "sample",
    col_names_WL
  )
  
  # skip columns if n_col is given
  if (!is.null(n_col)){
    
    range <- cell_limits(c(n_skip + 1, 1),
                          c(NA, n_col))
    data <- read_xlsx(path = file, col_names = col_names, range = range) 
  } else {
    
    data <- read_xlsx(path = file, col_names = col_names, skip = n_skip) 
    
  }
  
  
  
  # transform into long format and add new variables
  data <- data %>%
    gather(key = condition,
           value = f_signal,
           starts_with("WL")) %>%
    separate(condition, into = c("wave_length", "time"), remove = FALSE) %>%
    filter(wave_length %in%  str_c("WL", keep_wavelength)) %>% # filter only used wavelength
    mutate(timepoint = as.numeric(time),
           well_row = str_extract(well, "[A-H]"),
           well_col = str_extract(well, "\\d{2}"),
           experiment = exp,
           cell_type = cell_type,
           cell_number = cell_number,
           plate = as.character(plate),
           well_id = paste(experiment, plate, well, sep = "_")) %>% # added a variable that identify unique experimental reading
    select(-c(condition, time, sample, wave_length))
  
  
  if (plot_data){
    plot <- quick_plot(data)
    print(plot)
  }
  
  if (length(unique(data$well)) != 96){ # Raise a warning is less than 96 well are parsed (eg. if n_skip argument is incorrect)
    warning("number of parsed well is not 96")
  }
  
  if (!is.null(annotate)){
    if (annotate %in% c("minimal","wide")){
      data <- annotate_wells(df = data, type = annotate )
    }
  }
  return(data)
}

# Parse LDH results

parse_LDH <- function(file,
                      range = cell_limits(ul = c(14, NA), lr = c(NA, 4)),
                      col_names = c("well_row", "well_col", "sample", "LDH_assay"),
                      exp,
                      plate){
  
  
  
  # arguments for parsing
  dir  <- getwd()
  path <- str_c(
    dir,
    file
  )
  range <- range
  col_names <- col_names
  
  # new variables
  exp <- exp
  plate <- plate
  
  # print first 15 lines (for checking)
  out <- read_xlsx(file, n_max = 15)
  print(out)
  
  # parse data
  data <- read_xlsx(file, range = range, col_names = col_names) %>%
    mutate(well = str_c(well_row, formatC(well_col, width = 2, format = "d", flag = "0"))) %>%
    # discard unused variables
    select(-c(well_row, well_col, sample)) %>%
    mutate(experiment = exp, plate = plate)
  
  return(data)
  
  
}

# Parse growth curves

parse_growth <- function(file,
                         n_skip = 14,
                         well_cols = c("separated", "merged"), # are well columns separated in row-col (the default) or merged in one column?
                         interval = 15,
                         duration,
                         normalised = FALSE,
                         multiple_tranformations = NULL, # integer vector of length: first integer is the number of data transformations; second integer: the data transformation to choose from (e.g. raw data and blank normalised data - chose raw data: c(2,1))
                         exp,
                         plate,
                         plot_data = TRUE, # quick plot of raw data
                         annotate = c("minimal", "wide", "none")){
  
  annotate <- match.arg(annotate)
  well_cols <- match.arg(well_cols)
  
  # Starting messages
  print(glue::glue("Parsing growth data from file {file}"))
  
  
  # Timepoint calculation in minutes
  duration_h <- duration[1] * 60 # hour component of duration, transformed in minutes
  duration_min <- duration[2] # minute component of duration
  timepoint_max <- duration_h + duration_min
  # generates vector of timepoints (in minutes) of fluorescence measurement. Here: every 6 minutes until 22 hours and 18 minutes
  timepoints <- seq(from = 0, to = timepoint_max, by = interval) 
  
  
  
  # generate vector of column names based on timepoints 
  col_names_OD <- str_c("OD_", timepoints)
  
  if (!is.null(multiple_tranformations)){
    n_transf <- 1:multiple_tranformations[1]
    # print message 
    message(str_c("Raw data are provided in ", multiple_tranformations[1], " forms"))
    ordinal <- case_when(multiple_tranformations[2] == 1 ~ "st",
                         multiple_tranformations[2] == 2 ~ "nd",
                         multiple_tranformations[2] == 3 ~ "rd",
                         TRUE ~ "th")
    ordinal <- str_c(multiple_tranformations[2], ordinal)
    message(str_c("Only the ", ordinal, " form will be parsed" ))
    col_names_OD <- sapply(n_transf, function(x){
      str_c(x, "_", col_names_OD)
    })
  }
  
  if (well_cols == "merged") {
    col_names <- c("well","sample",col_names_OD)
  } else if (well_cols == "separated") {
    col_names <- c("well_row", "well_col", "sample",col_names_OD)
  }
  
  # parse data
  data <- read_xlsx(file, skip = n_skip, col_names = col_names) 
  
  # check that OD columns are numeric
  cols <- data %>%
    select(contains("OD_") & where(is.character)) %>%
    colnames()
  
  if (length(cols) > 0) {
    
    
    warning(glue::glue("{length(cols)} OD columns in wide format contain non-numeric values"))
    
    data <- data %>%
      mutate(across(all_of(cols), ~(if_else(is.na(as.numeric(.)),
                                            NA_real_,
                                            as.numeric(.)))))
    
  }
  
  # pivot in long form
  data <- data %>%
    pivot_longer(cols = contains("OD"), 
                 names_to = "time", 
                 values_to = "OD")
  
  if (!"well" %in% colnames(data)){
    data <- data %>%
      mutate(well = str_c(well_row, 
                          formatC(well_col, 
                                  width = 2, 
                                  format = "d", 
                                  flag = "0")))
  }
  
  if (!"well_row" %in% colnames(data)){
    data <- data %>%
      separate(col = well, into = c("well_row", "well_col"), sep = 1, remove = FALSE) %>%
      mutate(well_col = as.integer(well_col))
  }
  
  if (!is.null(multiple_tranformations)) {
    s <- str_c(multiple_tranformations[2], "_")
    data  <- data %>%
      filter(str_detect(time, s)) %>%
      mutate(time = str_remove(time, "\\d_"))
  }
  
  data <- data %>%
    mutate(time = str_remove(time, "OD_"),
           OD_blank_normalised = normalised,
           experiment = exp,
           plate = as.character(plate),
           well_id = paste(experiment, plate, well, sep = "_")) %>% # added a variable that identify unique experimental reading) 
    mutate(timepoint = as.numeric(time)) %>%
    select(-c(sample, time)) 
  
  
  if (plot_data){
    df_plot <- rename(data, f_signal = OD)
    plot <- quick_plot(df_plot)
    print(plot)
  }
  
  if (length(unique(data$well)) != 96){ # Raise a warning is less than 96 well are parsed (eg. if n_skip argument is incorrect)
    warning("number of parsed well is not 96")
  }
  
  if (!is.null(annotate)){
    if (annotate %in% c("minimal","wide")){
      data <- annotate_wells(df = data, type = annotate )
    }
  }
  return(data)
}

# Plot growth curves (this is temporary, we will write a general plotting function later on)




# Standardisation of data ----

# Cut experimental max duration to the shortest experiment 
cut_exp_duration <- function(df){
  shortest_exp_time <- df %>% 
    group_by(experiment_id) %>% 
    filter(timepoint == max(timepoint)) %>%
    ungroup() %>%
    filter(timepoint == min(timepoint)) %>%
    .$timepoint %>%
    unique(.)
  df <- df %>%
    filter(timepoint <= as.numeric(as.character(shortest_exp_time)))
  message(paste("Cutting to shortest experiment: ",  
              shortest_exp_time, 
              " min / ", 
              shortest_exp_time/60, 
              " hours"))
  return(df)
}

# Add a synchronised timepoint variable ($timepoint_sync) to kinetic data using JE2 point of max death rate as reference timepoint  
# require to run fit_smooth_spline first to get 1st derivative

sync_experiments <- function(df){
  
  if(!("fitted.deriv1" %in% colnames(df))) {
    stop("You need to run fit_smooth_spline on the data first")
  }
  
  # get point of max death rate for each well_id
  t <- df %>%
    group_by(well_id) %>%
    filter(fitted.deriv1 == max(fitted.deriv1)) %>%
    ungroup() %>%
    mutate(point_of_max_death_rate = timepoint, max_death_rate = fitted.deriv1) %>%
    select(well_id, point_of_max_death_rate, max_death_rate) %>%
    merge(df, ., by = "well_id") 
  
  # add mean of JE2 max death rate for each experiment_id to df for synchronisation
  t1 <- t %>%
    filter(sample_id == "JE2") %>%
    group_by(experiment_id) %>%
    summarise(JE2_point_of_max_death_rate = plyr::round_any(mean(point_of_max_death_rate), 6)) %>% # take mean of JE2 and round to a multipe of 6 (min)
    ungroup() %>%
    select(experiment_id, JE2_point_of_max_death_rate) %>%
    merge(t, ., by = "experiment_id")
  
  # calculate and apply JE2 time offset for each experiment_id using JE2 point of max death rate
  t2 <- t1 %>%
    mutate(time_offset = JE2_point_of_max_death_rate - min(JE2_point_of_max_death_rate)) %>%
    mutate(timepoint_sync = timepoint - time_offset) %>%
    merge(., t, by = "well_id") # FIXME
  
  return(t2)
  
}


# Deprecated standardisation function (replaced by strandardise_data() which standardise using fitted Je2 curve by default)
JE2_POMS_standardisation <- function(df, cut = T, plot = F){
  warning("This function is deprecated, you should use standardise_data() instead")
  if (cut){
    df <- cut_exp_duration(df)
  }
  
  print("Please check the number of JE2 replicate(s) per experiment_id that will be used for standardisation:")
  JE2_count_df <- df %>% 
    group_by(experiment_id, sample_id, replicate) %>% 
    filter(row_number() == 1) %>% 
    group_by(experiment_id, sample_id, .drop = F) %>%
    count(.drop = F) %>%
    filter(sample_id == "JE2") 
  print.data.frame(JE2_count_df)
  
  # Check if number of JE2 replicates in each experiment_id, 
  # if 0 replicate => experiment_id is excluded
  # if 1 replicate => raise a warning
  zero_rep_JE2_experiments <- JE2_count_df %>% filter(n == 0) %>% .$experiment_id
  one_rep_JE2_experiments <- JE2_count_df %>% filter(n == 1) %>% .$experiment_id
  if (length(zero_rep_JE2_experiments) != 0){
    warning("The following experiment(s) don't have any JE2 replicate and will be excluded: ", noBreaks. = F)
    warning(paste(zero_rep_JE2_experiments), noBreaks. = F)
    df <- df %>% filter(!(experiment_id %in% zero_rep_JE2_experiments))
  }
  if (length(one_rep_JE2_experiments) != 0){
    warning("The following experiment(s) have only one JE2 replicate that will be used for standardisation:")
    warning(paste(one_rep_JE2_experiments))
  }
  
  std_df <- df %>%
    filter(sample_id == "JE2") %>%
    group_by(experiment_id, timepoint) %>%
    summarise(f_signal_JE2.mean = mean(f_signal)) %>%
    ungroup() %>%
    group_by(experiment_id) %>% 
    summarise(f_signal_JE2.max = max(f_signal_JE2.mean),
              f_signal_JE2.min = min(f_signal_JE2.mean)) %>%
    ungroup() %>%
    merge(df, . , by = "experiment_id") %>%
    mutate(f_signal.std = (f_signal - f_signal_JE2.min) / (f_signal_JE2.max - f_signal_JE2.min))
  
  if (plot) {
    p <- plot_cell_death(std_df, standardized = T)
    print(p)
  }
  return(std_df)
}


# Perform proportion of maximum scaling standardisation with JE2
# Same as JE2_POMS_standardisation but now use mean of fitted JE2 signal to extract min and max
# This is to avoid standardising using outlier max an mean point from raw signal for standardisation (eg. machine reading error)
standardise_curves <- function(df, cut = T, plot = F){ 
  
  if (cut){
    df <- cut_exp_duration(df)
  }
  
  print("Please check the number of JE2 replicate(s) per experiment that will be used for standardisation:")
  # JE2_count_df <- df %>% 
  #    mutate(plate = as.character(plate)) %>%
  #   group_by(experiment, plate, sample_id, replicate) %>% 
  #   filter(row_number() == 1) %>% 
  #   droplevels() %>%
  #   group_by(experiment, plate, sample_id, .drop = F) %>%
  #   count(.drop = F) %>%
  #   filter(sample_id == "JE2") 
  JE2_count_df <- df %>%
    filter(sample_id == "JE2") %>%
    select(experiment_id, plate, sample_id, well_id) %>%
    distinct() %>%
    group_by(experiment_id, plate) %>%
    mutate(n = n()) %>%
    select(-well_id) %>%
    distinct() %>%
    ungroup() 
  print.data.frame(JE2_count_df)
  
  # Check if number of JE2 replicates in each experiment_id, 
  # if 0 replicate => experiment_id is excluded
  # if 1 replicate => raise a warning
  zero_rep_JE2_experiments <- JE2_count_df %>% filter(n == 0) %>% .$experiment_id %>% as.character()
  one_rep_JE2_experiments <- JE2_count_df %>% filter(n == 1) %>% .$experiment_id %>% as.character()
  if (length(zero_rep_JE2_experiments) != 0){
    warning("The following experiment(s) don't have any JE2 replicate and will be excluded: ", noBreaks. = F)
    warning(paste(c("\n", paste0(zero_rep_JE2_experiments, sep = "\n"), "\n")), noBreaks. = F)
    df <- df %>% filter(!(experiment_id %in% zero_rep_JE2_experiments))
  }
  if (length(one_rep_JE2_experiments) != 0){
    warning("The following experiment(s) have only one JE2 replicate that will be used for standardisation: ", noBreaks. = F)
    warning(paste(c("\n", paste0(one_rep_JE2_experiments, sep = "\n"), "\n")), noBreaks. = F)
  }
  
  std_df <- df %>%
    filter(sample_id == "JE2") %>%
    droplevels() %>%
    fit_curves(use_standardised = F, message = F) %>%                             # calculate fitted value for JE2
    group_by(experiment_id, timepoint) %>%
    summarise(f_signal_JE2.mean = mean(f_signal.fitted)) %>%                      # use fitted data to calculate means
    ungroup() %>%
    group_by(experiment_id) %>% 
    summarise(f_signal_JE2.max = max(f_signal_JE2.mean),
              f_signal_JE2.min = min(f_signal_JE2.mean)) %>%
    ungroup() %>%
    merge(df, . , by = "experiment_id") %>%
    mutate(f_signal.std = (f_signal - f_signal_JE2.min) / (f_signal_JE2.max - f_signal_JE2.min))
  
  if (plot) {
    p <- plot_cell_death(std_df, standardized = T)
    print(p)
  }
  return(std_df)
}

# function to filter/reduce number of point plotted when all_replicates = T
divisible <- function(x, y) {
  res = c()
  for (val in x) {
    if (val %% y == 0) {
      res = c(res, T)
    }else{
      res = c(res, F)
    }
  }
  return(res)
}

# Outliers detection functions (inspired from https://www.r-bloggers.com/combined-outlier-detection-with-dplyr-and-ruler/) ----

# different method to test if a value is an outlier

# Z-score outliers: outlier if the absolute difference with the mean is > to k times standard deviation
is_out_sd <- function(x, k_sd = 3, min_sd = 0, na.rm = TRUE) {
  (abs(x - mean(x, na.rm = na.rm)) > k_sd * sd(x, na.rm = na.rm)) & (sd(x) > min_sd)
}

# Z-score with median and Median Absolute Deviation (MAD):  outlier if absolute difference with the median is > to k times median absolute deviation
is_out_mad <- function(x, k_mad = 3, min_sd = 0, na.rm = TRUE) {
  (abs(x - median(x, na.rm = na.rm)) > k_mad * mad(x, na.rm = na.rm)) & (sd(x) > min_sd)
}

# Tukey's fences method = box plot outlier detection method based on interquartile range: outlier if x distance is > k times (default 1.5) 1st and last quartile   
is_out_tukey <- function(x, k_tukey = 1.5, min_sd = 0,  na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  ((x < quar[1] - k_tukey * iqr) | (x > quar[2] + k_tukey * iqr)) & (sd(x) > min_sd)
}

# Function to add outliers TRUE/FALSE variable for each time point in kinetic result dataframe (using the different method above)
flag_outlier_timepoints <- function(df, k_sd = 3, k_mad = 3, k_tukey = 1.5, min_sd = 0.2, use_standardised = T) { 
  
  if (use_standardised == F) {df <- df %>% mutate(f_signal.std = f_signal)} # used for growth curve data
  
  warning("Flagging outlier timepoints", "\n", noBreaks. = F)
  df <- df %>% 
    group_by(sample_id, timepoint) %>%
    mutate(sd = sd(f_signal.std)) %>%
    mutate(outlier_sd = is_out_sd(f_signal.std, k_sd = k_sd, min_sd = min_sd)) %>% # not sure if raw or standardised or fitted value is the most appropriate here
    mutate(outlier_mad = is_out_mad(f_signal.std, k_mad = k_mad, min_sd = min_sd )) %>% 
    mutate(outlier_tukey = is_out_tukey(f_signal.std, k_tukey = k_tukey, min_sd = min_sd)) %>%
    mutate(sd = ifelse(is.na(sd), 0, sd)) %>%
    mutate(outlier_sd = ifelse(is.na(outlier_sd), FALSE, outlier_sd)) %>%
    mutate(outlier_tukey = ifelse(is.na(outlier_tukey), FALSE, outlier_tukey)) %>%
    ungroup()
  
  if (use_standardised == F) {df <- df %>% mutate(f_signal = f_signal.std) %>% select(-f_signal.std)}
  
  return(df)
}

summarise_outlier <- function(df, 
                              k_sd = 3,          # Set thresholds for outlier detection 
                              k_mad = 3, 
                              k_tukey = 1.5,
                              min_sd = 0.2) 
                             
{
  # flag outliers if not previously done
  if(!("outlier_sd" %in% colnames(df))) {
    df <- flag_outlier_timepoints(df, k_sd = k_sd, k_mad = k_mad, k_tukey = k_tukey, min_sd = min_sd )
  }else{
    warning("Skipping outlier flagging with flag_outlier_timepoints() and use previously flagged outlier timepoint", noBreaks. = F)
  }
  outlier_df <- df %>%
    group_by(well_id)%>%
    summarise(proportion_outlier_sd = sum(outlier_sd)/ length(outlier_sd),
              proportion_outlier_mad = sum(outlier_mad)/ length(outlier_mad),
              proportion_outlier_tukey = sum(outlier_tukey)/ length(outlier_tukey))
  
    # merge results with well metadata
    merge(df, outlier_df, by = "well_id") %>%
      merge(., df %>%
              filter(!duplicated(well_id)) %>% 
              select(-starts_with("f_signal"), 
                     -starts_with("timepoint"), 
                     -starts_with("fitted"), 
                     -starts_with("outlier")),
           by = "well_id", 
           all.x = T, 
           all.y = F) 
}

flag_outlier_curves <- function(df, 
                                use_standardised = T,
                                k_sd = 3,                         # Set thresholds for outlier detection 
                                k_mad = 3, 
                                k_tukey = 1.5,
                                min_sd = 0.2,
                                outlier_sd_proportion = 0.1,      # Set the minimum proportion of outlier timepoint in a curve (well_id) to consider outlier
                                outlier_mad_proportion = 0.1,
                                outlier_tukey_proportion = 0.1,
                                plot = F,                         # Set to True to visualise outliers point curves
                                filter = F)                       # Set to True to remove flagged outliers from the dataset
  
{
  # flag outliers if not previously done
  if(!("outlier_sd" %in% colnames(df))) {
    df <- flag_outlier_timepoints(df, k_sd = k_sd, k_mad = k_mad, k_tukey = k_tukey, min_sd = min_sd, use_standardised = use_standardised)
  }else{
    warning("Skipping outlier flagging with flag_outlier_timepoints() and use previously flagged outlier timepoint", noBreaks. = F)
  }
  
  # calculate proportion of outlier timepoint for each curve
  df_out <- df %>%
    group_by(well_id)%>%
    mutate(proportion_outlier_sd = sum(outlier_sd)/ length(outlier_sd),
           proportion_outlier_mad = sum(outlier_mad)/ length(outlier_mad),
           proportion_outlier_tukey = sum(outlier_tukey)/ length(outlier_tukey))
  

  # label outlier curves according to proportion threshold
  df_out <- df_out %>% mutate(outlier_curve = proportion_outlier_sd > outlier_sd_proportion |
                                proportion_outlier_mad > outlier_mad_proportion |
                                proportion_outlier_tukey > outlier_tukey_proportion)
  
  # df with filtered outliers
  df_filt <- df_out %>% filter(outlier_curve == F)
  
  # count well and outliers
  nb_well = length(unique(df_out$well_id))
  nb_outlier = nb_well - length(unique(df_filt$well_id))
  percent_outlier = (nb_outlier/nb_well) * 100
  
  outlier_curve.df <- df_out %>%
    filter(outlier_curve == T) %>%
    filter(!duplicated(well_id)) %>%
    select(well_id, sample_id, experiment_id, plate) %>%
    arrange(well_id)
  
  
  # print counts
  
  print(paste("Total number of curves:", as.character(nb_well)))
  print(paste("Total number of outliers detected:", as.character(nb_outlier)))
  print(paste("Percentage of outliers detected:", as.character(percent_outlier), "%"))
  print("Outliers curves:")
  print.data.frame(outlier_curve.df)
  
  # plot data
  if (plot) {

    plot_cell_death(df = df_out,
                    standardized = T,
                    fitted = T,
                    combine_experiments = "all", 
                    comparator = "JE2",
                    highlight_well1 = outlier_curve.df$well_id,
                    highlight_point1 = "outlier_mad",
                    highlight_point2 = "outlier_tukey",
                    highlight_point3 = "outlier_sd")
  }
  
  # return filtered or unfiltered data
  if (filter){
    return(df_filt)
  }else{
    return(df_out)
  }
}




# Function to check JE2
check_reference <- function(df,
                            reference = "JE2",                # Set sample_id of reference
                            outlier_sd_proportion = 0.05,      # Set the proportion of oultier timepoint in a curve (well_id) to consider outlier
                            outlier_mad_proportion = 0.05,
                            outlier_tukey_proportion = 0.05,
                            k_sd = 3,                         # Set thresholds for outlier detection at each timepoint
                            k_mad = 3, 
                            k_tukey = 1.5,
                            min_sd = 0.2,                     # Set min sd of sample to consider to retain oulier timepoint
                            plot = T,
                            filter = F)                     
{
 message("Identifying ", reference, " outliers:")
 df_ref <- df %>% 
   filter(sample_id == reference) %>%
   droplevels() %>%
   standardise_curves() %>%
   fit_curves() %>%
   flag_outlier_curves(k_sd = k_sd,                                        # Set thresholds for outlier detection 
                       k_mad = k_mad, 
                       k_tukey = k_tukey,
                       min_sd = min_sd,
                       outlier_sd_proportion = outlier_sd_proportion,      # Set the minimum proportion of outlier timepoint in a curve (well_id) to consider outlier
                       outlier_mad_proportion = outlier_mad_proportion,
                       outlier_tukey_proportion = outlier_tukey_proportion,
                       plot = plot,                                        # Set to True to visualise outliers point curves
                       filter = filter) %>%
   ungroup()
}  
                            

# Model fitting ----

# fit smooth spline to each PI measurement (well_id)

fit_smooth_spline <- function(df, spar = 0.8, growth_curves = F){
  warning("This function is deprecated, you should use fit_curves()")
  df <- df %>% droplevels() # Avoid subscript out of bounds error when data have been filtered with a factor variable (eg. experiment)
  
  if (growth_curves) {
    df <- df %>% rename(f_signal = OD)
    message("Fitting OD signal for growth curve")
  } else {
    message("\n","Fitting PI signal for cell death")
  }
  
  df_res <- df
  
  if(!("f_signal.std" %in% colnames(df))) {
    warning("Using non-standardised signal for fitting", noBreaks. = F)
    warning("Run JE2_POMS_standardisation() function first if you want to use standardised signal", noBreaks. = F)
    df_fit <- df %>% 
      group_by(well_id) %>% 
      do(fit = augment(smooth.spline(x = .$timepoint, y = .$f_signal, spar = spar, keep.data = T))) %>%
      unnest(cols = c(fit)) %>%
      ungroup() %>%
      transmute(well_id = well_id, timepoint = x, f_signal = y, f_signal.fitted = .fitted, f_signal.residual = .resid) %>%
      select(-f_signal) %>%
      merge(., df_res, by = c("well_id", "timepoint"))
    
    # calculate 1st and 2nd derivatives
    message("Calculating 1st and 2nd derivatives")
    df_fit <- df_fit %>%
      group_by(well_id) %>%
      do(timepoint = predict(smooth.spline(x = .$timepoint, y = .$f_signal, spar = spar, keep.data = T), deriv = 0)$x,
         fitted.deriv1 = predict(smooth.spline(x = .$timepoint, y = .$f_signal, spar = spar, keep.data = T), deriv = 1)$y,
         fitted.deriv2 = predict(smooth.spline(x = .$timepoint, y = .$f_signal, spar = spar, keep.data = T), deriv = 2)$y) %>%
      unnest(cols = c(timepoint, fitted.deriv1, fitted.deriv2 )) %>%
      ungroup() %>%
      merge(df_fit, .,  by = c("well_id", "timepoint"))
    
  }else{
    message("Using standardised signal for fitting")
    df_fit <- df %>% 
      group_by(well_id) %>% 
      do(fit = augment(smooth.spline(x = .$timepoint, y = .$f_signal.std, spar = spar, keep.data = T))) %>% 
      # run smooth.spline and store the resulting list in a variable fit
      #broom::augment(fit) %>% # add components of the list as new variables
      unnest(cols = c(fit)) %>%
      ungroup() %>%
      transmute(well_id = well_id, timepoint = x, f_signal.std = y, f_signal.fitted = .fitted, f_signal.residual = .resid) %>%
      select(-f_signal.std) %>%
      merge(., df_res, by = c("well_id", "timepoint"))
    
    # calculate 1st and 2nd derivatives
    message("Calculating 1st and 2nd derivatives")
    df_fit <- df_fit %>%
      group_by(well_id) %>%
      do(timepoint = predict(smooth.spline(x = .$timepoint, y = .$f_signal.std, spar = spar, keep.data = T), deriv = 0)$x,
         fitted.deriv1 = predict(smooth.spline(x = .$timepoint, y = .$f_signal.std, spar = spar, keep.data = T), deriv = 1)$y,
         fitted.deriv2 = predict(smooth.spline(x = .$timepoint, y = .$f_signal.std, spar = spar, keep.data = T), deriv = 2)$y) %>%
      unnest(cols = c(timepoint, fitted.deriv1, fitted.deriv2 )) %>%
      ungroup() %>%
      merge(df_fit, .,  by = c("well_id", "timepoint"))
  }
  
  if (growth_curves) df_fit <- rename(df_fit, OD = f_signal, OD.fitted = f_signal.fitted, OD.residual = f_signal.residual)
  return(df_fit)

}

# fit smooth spline to each PI measurement (well_id)
# !!! This new function doesn't calculate 1st and 2nd drivative anymore --> this is done by get_parameters() now
fit_curves <- function(df, spar = 0.8, use_standardised = T, growth_curves = F, message = T){
  df <- df %>% droplevels() # Avoid subscript out of bounds error when data have been filtered with a factor variable (eg. experiment)
  
  if (growth_curves) {
    df <- df %>% rename(f_signal = OD)
    if (message) message("\n","Fitting OD signal for growth curve")
  } else {
    if (message) message("\n","Fitting PI signal for cell death", appendLF = F)
  }
  
  df_res <- df
  
  if(use_standardised){
    if (message) message(" using standardised data")
    df_fit <- df %>% 
      group_by(well_id) %>% 
      do(fit = augment(smooth.spline(x = .$timepoint, y = .$f_signal.std, spar = spar, keep.data = T))) %>% 
      # run smooth.spline and store the resulting list in a variable fit
      # broom::augment(fit) %>% # add components of the list as new variables
      unnest(cols = c(fit)) %>%
      ungroup() %>%
      transmute(well_id = well_id, timepoint = x, f_signal.std = y, f_signal.fitted = .fitted, f_signal.residual = .resid) %>%
      select(-f_signal.std) %>%
      merge(., df_res, by = c("well_id", "timepoint"))
    
  }else{
    
    if (message) message(" using non-standardised data")
    #warning("Run standardise_curves() function first if you want to use standardised signal", noBreaks. = F)
    df_fit <- df %>% 
    group_by(well_id) %>% 
    do(fit = augment(smooth.spline(x = .$timepoint, y = .$f_signal, spar = spar, keep.data = T))) %>%
    #broom::augment(fit) %>%
    unnest(cols = c(fit)) %>%
    ungroup() %>%
    transmute(well_id = well_id, timepoint = x, f_signal = y, f_signal.fitted = .fitted, f_signal.residual = .resid) %>%
    select(-f_signal) %>%
    merge(., df_res, by = c("well_id", "timepoint"))
    
  }
  if (growth_curves) df_fit <- rename(df_fit, OD = f_signal, OD.fitted = f_signal.fitted, OD.residual = f_signal.residual)
  return(df_fit)
}

# Parameters extraction from curves ----

# get max death + timepoint of max death from smooth spline
get_point_max_death <- function(df){
  
  # if(!("f_signal.fitted" %in% colnames(df))) {
  #   stop("You need to run fit_smooth_spline on the data first")
  # }
  
  t <- df %>%
    group_by(well_id) %>%
    filter(f_signal.fitted == max(f_signal.fitted)) %>%
    ungroup() %>%
    mutate(time_of_max_death = timepoint, max_death = f_signal.fitted) %>%
    select(well_id, time_of_max_death, max_death)  
  return(t)
}

# get min death + timepoint of min death from smooth spline
get_point_min_death <- function(df){
  
  # if(!("f_signal.fitted" %in% colnames(df))) {
  #   stop("You need to run fit_curves on the data first")
  # }
  
  t <- df %>%
    group_by(well_id) %>%
    filter(f_signal.fitted == min(f_signal.fitted)) %>%
    ungroup() %>%
    mutate(time_of_min_death = timepoint, min_death = f_signal.fitted) %>%
    select(well_id, time_of_min_death, min_death)  
  return(t)
}

# get max death rate + timepoint of max death rate from smooth spline derivative (need to run fit_smooth spline first)
get_point_max_death_rate <- function(df){
  
   # if(!("fitted.deriv1" %in% colnames(df))) {
   #   stop("You need to run get_derivatives() on the data first to get 1st derivative (and death/growth rate)")
   # }
  
  t <- df %>%
    group_by(well_id) %>%
    filter(fitted.deriv1 == max(fitted.deriv1)) %>%
    ungroup() %>%
    mutate(time_of_max_rate_death = timepoint, max_rate_death = fitted.deriv1, doubling_time_death = log(2)/max_rate_death) %>%
    select(well_id, time_of_max_rate_death, max_rate_death, doubling_time_death)  
  return(t)
}

# get min death rate + timepoint of min death rate from smooth spline derivative (need to run fit_smooth spline first)
get_point_min_death_rate <- function(df){
  
  # if(!("fitted.deriv1" %in% colnames(df))) {
  #   stop("You need to run fit_curves on the data first to get 1st derivative (and death/growth rate)")
  # }
  
  t <- df %>%
    group_by(well_id) %>%
    filter(fitted.deriv1 == min(fitted.deriv1)) %>%
    ungroup() %>%
    mutate(time_of_min_rate_death = timepoint, min_rate_death = fitted.deriv1) %>%
    select(well_id, time_of_min_rate_death, min_rate_death)  
  return(t)
}

# get end point signal
get_end_point <- function(df){
  
  # if(!("fitted.deriv1" %in% colnames(df))) {
  #   stop("You need to run fit_curves on the data first to get 1st derivative (and death/growth rate)")
  # }
  
  t <- df %>%
    group_by(well_id) %>%
    filter(timepoint == max(timepoint)) %>%
    ungroup() %>%
    mutate(end_point_death = f_signal.fitted, end_point_timepoint_death = timepoint) %>%
    select(well_id, end_point_death, end_point_timepoint_death)
  return(t)
}

# get AUC (need to run fit_smooth spline first)
get_AUC <- function(df){
  
  # if(!("fitted.deriv1" %in% colnames(df))) {
  #   stop("You need to run fit_curves on the data first")
  # }
  
  # Calculating area under the curve (this quick way is is equivalent to grofit method and avoid fitting smooth splines twice)
  t <- df %>%
    group_by(well_id) %>%
    summarise(AUC_death = sum(f_signal.fitted)) %>%
    ungroup() %>%
    select(well_id, AUC_death)
  return(t)
}

get_derivatives <- function(df, spar = 0.8, use_fitted = F){
  
  used_data <- "standardised"
  
  if (use_fitted == T) {
    df <- df %>% mutate(f_signal.std = f_signal.fitted)
    used_data <- "fitted"
  }
  
  # calculate 1st and 2nd derivatives
  message(paste0(c("Calculating 1st and 2nd derivatives from ", used_data, " data")))
  df <- df %>%
    group_by(well_id) %>%
    do(timepoint = predict(smooth.spline(x = .$timepoint, y = .$f_signal.std, spar = spar, keep.data = T), deriv = 0)$x,
       fitted.deriv1 = predict(smooth.spline(x = .$timepoint, y = .$f_signal.std, spar = spar, keep.data = T), deriv = 1)$y,
       fitted.deriv2 = predict(smooth.spline(x = .$timepoint, y = .$f_signal.std, spar = spar, keep.data = T), deriv = 2)$y) %>%
    unnest(cols = c(timepoint, fitted.deriv1, fitted.deriv2 )) %>%
    ungroup() %>%
    merge(df, .,  by = c("well_id", "timepoint"))
}
  

# fit growth/death parameters
get_parameters <- function(df, growth_curves = F, spar = 0.8, use_fitted = F){
  
  
  if (growth_curves) df <- rename(df, f_signal.fitted = OD.fitted)
  
  if(!("f_signal.fitted" %in% colnames(df))) {
    warning("fit_curves() is running\n")
    df <- fit_curves(df, spar = spar)
  }
  
  if(!("fitted.deriv1" %in% colnames(df))) {
     warning("get_derivatives() is running\n")
     df <- get_derivatives(df, spar = spar, use_fitted = use_fitted)
  }
  
  if (use_fitted == T) df <- df %>% mutate(f_signal.std = f_signal.fitted) # used for growth curve data
  
  
  # results with results from get_point_max_death_rate() 
  res.df <- get_point_max_death_rate(df)
  
  # merge results from get_AUC() 
  res.df <- merge(res.df,
                  get_AUC(df),
                  by = "well_id", all.x = T, all.y = F)
  
  # merge results with results from get_point_max_death() 
  res.df <- merge(res.df,
                  get_point_max_death(df),
                  by = "well_id", all.x = T, all.y = F)
  
  # merge results with results from get_point_min_death() 
  res.df <- merge(res.df,
                  get_point_min_death(df),
                  by = "well_id", all.x = T, all.y = F)
  
 
  
  # merge results with results from get_point_min_death_rate() 
  #res.df <- merge(res.df,
  #                get_point_min_death_rate(df),
  #                by = "well_id", all.x = T, all.y = F)
  
  # merge results with results from get_end_point
  res.df <-  merge(res.df,
                   get_end_point(df),
                   by = "well_id", all.x = T, all.y = F)
                  
  # merge results with well metadata
  res.df <- merge(res.df, df %>% 
                    filter(!duplicated(well_id)) %>% 
                    select(-starts_with("f_signal"), -starts_with("timepoint"), -starts_with("fitted")),
                  by = "well_id", 
                  all.x = T, 
                  all.y = F) 
  return(res.df)
}

# Use ggpubr to test for significant differences in growth/death parameters
get_parameters_stat <- function(df_param,
                               ref = "JE2",
                               test = "wilcox.test",
                               group = NULL) { # a character vector containing the name of grouping variables (eg. c("experiment"))
  require(ggpubr)
  t <- rbind(compare_means(formula =  max_rate_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  time_of_max_rate_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             # compare_means(formula =  min_rate_death ~ sample_id, 
             #               data = df_param, 
             #               ref.group = ref, 
             #               method = test,
             #               group.by = group),
             # compare_means(formula =  time_of_min_rate_death ~ sample_id, 
             #               data = df_param, 
             #               ref.group = ref, 
             #               method = test,
             #               group.by = group),
             compare_means(formula =  max_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  time_of_max_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  min_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  time_of_min_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  AUC_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group))
  colnames(t)[1] <- "tested_parameter"
  colnames(t)[2] <- "reference"
  colnames(t)[3] <- "sample_id"
  
  if (!is.null(group)){
    colnames(t)[1] <- "grouping"
    colnames(t)[2] <- "tested_parameter"
    colnames(t)[3] <- "reference"
    colnames(t)[4] <- "sample_id"
    
  }
  return(t)
}

## Within-host evolution analysis functions
pair_comparator <- function(df_pairs, # a dataframe with standardised and fitted toxicity data + pair metadata
                            save_plot = F,
                            print_plot = T,
                            return_plot = F){
  type = "fitted smooth spline"
  pairs <- df_pairs %>% filter(strain_group != "CONTROL") %>% .$strain_group %>% unique()
  
  for (i in pairs) {
    my_plate <- df_pairs %>%
      filter(strain_group == i) %>%
      .$plate %>%
      unique()
    
    my_group <- c(i, "CONTROL")
    
    my_df <- df_pairs %>%
      filter(plate == my_plate) %>%
      filter(strain_group %in% my_group) %>%
      mutate(sample_id = if_else(
        strain_group == i,
        str_c(sample_id, "\n(", sample_type, ")"),
        sample_id
      ))
    
    my_comparator <- my_df %>%
      filter(intrahost_index == 1) %>%
      .$sample_id %>%
      unique()
    
    my_st <- unique(my_df %>% drop_na(ST) %>% .$ST %>% unique())
    my_pair <- str_c(i, " / ST ", my_st)
    
    my_title = str_c("Plate ", my_plate, " - ", my_pair, " - ", type)
    
    p <- plot_cell_death(df = my_df,
                         fitted = T,
                         combine_experiments = "plate", 
                         comparator = my_comparator,
                         all_replicates = 24,
                         title = my_title,
                         print_plot = FALSE,
                         return_plots = TRUE)
    
    my_metadata <- my_df %>%
      filter(intrahost_index == 0) %>%
      mutate(genes_short = str_wrap(genes, 50)) %>%
      # mutate(genes_short = if_else(str_length(genes) >= 50,
      #                           insert_new_line(genes, 50, check_length = FALSE),
      #                           genes)) %>%
      select(
        sample = sample_id,
        delay = intrahost_sampledelay,
        persister_type,
        scv,
        vancodiff = vanco_difference,
        mut = mut_count,
        genes = genes_short
      ) %>%
      distinct()
    
    t <- tableGrob(my_metadata)
    
    p <- p / t
    if (print_plot) print(p)
    
    if (save_plot != F){
      
      # check directory
      if (is.character(save_plot)){
        my_dir <- save_plot
      } else {
        my_dir <- "figures"
      }
      
      if (!dir.exists(dir)) dir.create(dir) # check if directory already exists (to avoid error)
      
      # create file name
      file <- str_c(dir, "/PI_kinetics_", i, "_", type, "_with_", unique(my_df$cell_number), "_cells_",
                    my_df$plate,".pdf")
      
      # calculate size based on number of graphs ( width = 3.5 in per graph)
      my_width <- n_distinct(my_df$sample_id) * 2
      
      ggsave(file, width = my_width, height = 5)
    }
    
  }
}

plot_parameters <- function(df, 
                            parameter, 
                            x = "sample_id",
                            label = "", 
                            order = NULL, 
                            refgroup = "JE2",
                            do_stat = T,  
                            method = "wilcox.test") {
  p <- ggboxplot(data = df,
                 y = parameter,
                 x = x,
                 color = "strain_group",
                 order = order ) +
    labs(x = "", y = label) +
    theme_bw() +
    theme(legend.position = "none")
  
  if (do_stat) {
    p = p +
      stat_compare_means(ref.group = refgroup,
                         method = method,
                         label = "p.signif")
  }
  
  return(p)
  
}

pair_comparator_parameters <- function(df_param, # a dataframe of cell death paramaters + pairs metadata
                                       parameters = "all", 
                                       method = "wilcox.test",
                                      save_plot = F,
                                       print_plot = T,
                                      return_plot = F){
  
  # fix parameters
  if (parameters == "all"){
    parameters <- c("AUC_death", "max_death", "max_rate_death", "time_of_max_rate_death")
    names(parameters) <- c("AUC\n",
                           "Peak cell death\n",
                           "Peak cell death rate (slope)\n",
                           "Timepoint of peak cell death rate\n")
  }
  
  pairs <- df_param %>% filter(strain_group != "CONTROL") %>% .$strain_group %>% unique()
  
  require(ggpubr)
  require(patchwork)
  for (i in pairs) {
    
    my_plate <- df_param %>%
      filter(strain_group == i) %>%
      .$plate %>%
      unique()
    
    my_group <- c(i, "CONTROL")
    
    my_df <- df_param %>%
      filter(plate == my_plate) %>%
      filter(strain_group %in% my_group) %>%
      mutate(sample_id = if_else(
        strain_group == i,
        str_c(sample_id, "\n(", sample_type, ")"),
        sample_id
      ))
    
    my_refgroup <- my_df %>%
      filter(intrahost_index == 1) %>%
      .$sample_id %>%
      unique()
    
    my_st <- unique(my_df %>% drop_na(ST) %>% .$ST %>% unique())
    
    my_pair <- str_c(i, " / ST ", my_st)
    
    my_title = str_c("Cell death parameters - ", my_pair, " (plate  ", my_plate, ")")
    my_subtitle = str_c("Comparator: ", str_remove(my_refgroup, "\\(index\\)"))
    
    my_order <- my_df %>%
      arrange(sample_id) %>%
      .$sample_id %>%
      unique()
    
    plots <- purrr::imap(parameters, plot_parameters, refgroup = my_refgroup, order = my_order, df = my_df,
                         method = method)
    
    p <- wrap_plots(plots) +
      plot_annotation(title = my_title,
                      subtitle = my_subtitle)
    
    if (print_plot) print(p)
    
    if (save_plot != F){
      
      # check directory
      if (is.character(save_plot)){
        my_dir <- save_plot
      } else {
        my_dir <- "figures"
      }
      
      if (!dir.exists(dir)) dir.create(dir) # check if directory already exists (to avoid error)
      
      file <- str_c(my_dir,
                    "/Cell_death_parameters_",
                    i,
                    "_with_",
                    unique(my_df$cell_number),
                    "_cells.pdf")
      
      ggsave(file, width = 9, height = 6)
    }
    
  }
}


# Plotting functions ----

# plot PI cell death function
plot_cell_death <- function(df, 
                            standardized = F, # set to TRUE to plot stadardised data (need to run JE2_POMS_standardisation first)
                            sync = F, # set to TRUE to synchronised timepoints according to JE2 max tox rate (need to run JE2_sync first)
                            fitted = F, # set to TRUE to plot fluorescence infered from smooth spline models (need to run fit_smooth_spline first)
                            residuals = F,  # set to TRUE to plot residuals from smooth spline models (need to run fit_smooth_spline first)
                            all_replicates = 24, # set to FALSE if you only want mean + standard deviation and not all replicate points/curves can be set to a numeric (eg. 24 will plot value every 24 min only)
                            combine_experiments = c("all"), # can set to any string corresponding to a column/variable, "all" to plot the combined results from different experiments, "plate" to combine results from different experiments for the same plate
                            annotate = F, # annotate well_id on plot (eg. to identify well_id of outliers)
                            strain_description = F, # add short description of the strain 
                            comparator = NULL, # add a comparator strain to each facet (set comparator to a sample_id)
                            highlight_well1 = c(), # set to a vector of well_id that you want to label/highlight
                            highlight_well2 = c(), # as above but with a different colour
                            highlight_point1 = NULL, # set to a string that correspond to a boolean (TRUE/FALSE) variable/column in your data (eg. "outlier_tukey" generated by flag_outlier_timepoints) 
                            highlight_point2 = NULL, # as above but with a different colour
                            highlight_point3 = NULL, # as above but with a different colour
                            derivative1 = F, # plot 1st derivatives
                            derivative2 = F, # plot 2nd derivatives
                            nrow = NULL, # number of plot rows passed to facet_wrap 
                            ncol = 8, # number of plot columns passed to facet_wrap 
                            ylim = NULL,
                            xlim = NULL,
                            print_plot = T, # set to FALSE to skip plot printing (can take time)
                            save_plot = F, # if not F specify directory to save the plot to
                            plot_size = NULL, # size (width, height) of the saved plot a numeric vector (if not supplied uses the size of the current graphics device)
                            return_plots = F, # return list of plots (or single plot) if set to TRUE
                            title = NULL, #set title to NULL to allow adding the title manually (default: string with plate number and data description)
                            subtitle = NULL) { #set subtitle (default: empty)
  

  
  # Check function parameters and set plot type
  combine_experiments <- match.arg(combine_experiments, choices = c("all", "plate", "experiment"), several.ok = F) # use match.arg() to extract the default ("none") if no value is given
  
  if (combine_experiments == "all") {
    df$comb <- "all"
    } else if (combine_experiments == "plate") {
    df$comb <- df$plate
    } else {
    df$comb <- df[,combine_experiments]
  }
  
  if (standardized) {
    df <- df %>%
      mutate(f_signal = f_signal.std)
    type = "standardised"
  }else{
    type = "non-standardised"
  }
  
  if (fitted) {
    df <- df %>%
      mutate(f_signal = f_signal.fitted)
    type = "fitted smooth spline"
  }
  
  if (residuals) {
    df <- df %>%
      mutate(f_signal = f_signal.residual)
    type = "residual of fitted smooth spline" 
  }
  
  if (derivative1) {
    df <- df %>%
      mutate(f_signal = fitted.deriv1)
    type = "fitted 1st derivative" 
  }
  
  if (derivative2) {
    df <- df %>%
      mutate(f_signal = fitted.deriv2)
    type = "fitted 2nd derivative" 
  }
  
  if (sync) {
    df <- df %>%
      mutate( timepoint = timepoint_sync) %>%
      cut_exp_duration() %>%
      filter(timepoint >= 0)
    type = str_c(type, " - synchronised")
  }

  # Main code to generate the plot
  require("Hmisc")
  all_p = list()
  for (exp in unique(df$comb)) {
    exp <- as.character(exp)
    
    df_exp <- df %>%
      filter(comb == exp)
    
    if (is.character(title)) {
      my_title <- title
    } else if (combine_experiments == "plate") {
      my_title <- str_c("Combined Plate ", exp,  " - ",  type)
    } else {
      my_plates <- str_c(unique(df_exp$plate), collapse = "+")
      my_title <- str_c(exp, " Plate ", my_plates, " - ",  type)
    }
    
    if (!is.null(subtitle)){
      my_subtitle <- subtitle
    } else {
      my_subtitle <- NULL
    }
    
    if (strain_description) {
      facets <- c("sample_id", "strain_description")
    } else {
      facets <- "sample_id"
    }
      
    p <- df_exp %>%
      ggplot(aes(x = timepoint/60, y = f_signal)) +
      #stat_summary(geom = "ribbon", fun.data = mean_cl_normal, fun.args=list(conf.int=0.95), fill = "#045a8d", alpha = 0.2) + # Manually checked and give same results, but 95% CI appear to be too narrow 
      stat_summary(geom = 'ribbon', fun.data = 'mean_sdl', fun.args = list(mult = 1), fill = "#045a8d", alpha = 0.2) + # Use standard deviation instead
      stat_summary(geom = "point", fun = "mean", color = "#045a8d", size = 1) +
      facet_wrap(facets, nrow = nrow, ncol = ncol) +
      labs(title = my_title,
           subtitle = my_subtitle,
           x = "Time (hour)",
           y = "Fluorescence intensity") +
      theme_bw()
    
    print(glue::glue("Basic plot for {my_title} created"))
    
    if (annotate) {
      require(ggrepel)
      p <- p + 
        geom_text_repel(data = subset(df_exp,
                                      timepoint == (count(df_exp, timepoint) %>% 
                                                      filter(n == max(n)) %>% 
                                                      filter(timepoint == max(timepoint)) %>% .$timepoint)),
                        aes(label = well_id),
                        force = 1,
                        max.iter = 10,
                        segment.size  = 0.5,
                        segment.color = "black",
                        direction     = "both")
    }
    
    if (!is.null(comparator)) {
      comparator_df <- df_exp %>%
        filter(sample_id == comparator) %>%
        select(-c(sample_id, strain_description))
      p <- p + 
        stat_summary(data = comparator_df, geom = 'ribbon', fun.data = 'mean_sdl', fun.args = list(mult = 1), fill = "#b30000", alpha = 0.2) +
        stat_summary(data = comparator_df, geom = "point", fun = "mean", color = "#b30000", alpha = 0.2, size = 1)
    }
    
    if (all_replicates != F) {
      if (is.numeric(all_replicates)) { # if all_replicate is a numeric set fitering value accordingly
        filt_div = all_replicates
      }else{
        filt_div = 6 # plot measure every 6 min by default (all_replicate = T)
      }
      p <- p +
        geom_point(data = df_exp %>% filter(divisible(timepoint, filt_div)), size = .05, color = "#045a8d")
      if (!is.null(comparator)) {
        p <- p + 
          geom_point(data = comparator_df %>% filter(divisible(timepoint, filt_div)), size = .05, color = "#b30000")
      }
    }
    
    if (!is_empty(highlight_well1)) {
      require(ggrepel)
      p <- p + 
        geom_text_repel(data = df_exp %>%
                          filter(well_id %in% highlight_well1) %>%
                          subset(., timepoint == (count(., timepoint) %>%
                                                    filter(n == max(n)) %>% 
                                                    filter(timepoint == max(timepoint)) %>% 
                                                    .$timepoint)),
                        aes(label = well_id),
                        max.iter = 10,
                        force = 4,
                        segment.size  = 0.5,
                        segment.color = "#4dac26",
                        colour = "#4dac26",
                        #size = 1,
                        direction = "both")
    }
    
    if (!is_empty(highlight_well2)) {
      require(ggrepel)
      p <- p + 
        geom_text_repel(data = df_exp %>%
                          filter(well_id %in% highlight_well2) %>%
                          subset(., timepoint == (count(., timepoint) %>%
                                                    filter(n == max(n)) %>% 
                                                    filter(timepoint == max(timepoint)) %>% 
                                                    .$timepoint)),
                        aes(label = well_id),
                        max.iter = 10,
                        force = 4,
                        segment.size  = 0.5,
                        segment.color = "#d01c8b",
                        colour = "#d01c8b",
                        #size = 1,
                        direction = "both")
    }
    
    if (is.character(highlight_point1)) {
      require(ggrepel)
      require(rlang)
      column_name <- rlang::sym(highlight_point1)
      p <- p + 
        geom_point(data = df_exp %>% filter(UQ(column_name) == T), colour = "#e6550d")
    }
    
    if (is.character(highlight_point2)) {
      require(ggrepel)
      require(rlang)
      column_name <- rlang::sym(highlight_point2)
      p <- p + 
        geom_point(data = df_exp %>% filter(UQ(column_name) == T), colour = "#d01c8b")
    }
    
    if (is.character(highlight_point3)) {
      require(ggrepel)
      require(rlang)
      column_name <- rlang::sym(highlight_point2)
      p <- p + 
        geom_point(data = df_exp %>% filter(UQ(column_name) == T), colour = "#4dac26")
      
    }
    
    # Set y and x limimts of plot
    if (!is_empty(ylim)) {
      p <- p + ylim(ylim)
    }
    
    if (!is_empty(xlim)) {
      p <- p + ylim(xlim)
    }
    
    # Save the plot
    
    if (save_plot != F) {
      
      # check directory
      if (is.character(save_plot)) {
        my_dir <- save_plot
      } else {
        my_dir <- "figures"
      }
      if (!dir.exists(my_dir))
        dir.create(my_dir) # check if directory already exists (to avoid error)
      
      # create file name
      my_title_file <- str_replace_all(my_title, "\\s+|-", "_") %>%
        str_replace_all("(?<=_)_", "")
      file <- str_c(
        my_dir,
        "/PI_kinetics_",
        my_title_file,
        "_with_",
        unique(df_exp$cell_number),
        "_cells.pdf"
      )
      
      # check if plot size is given as an argument
      if (is.null(plot_size)) {
        my_width <- NA
        my_height <- NA
      } else {
        my_width <- plot_size[1]
        my_height <- plot_size[2]
      }
      
      # save
      ggsave(file,
             width = my_width,
             height = my_height,
             limitsize = FALSE)
    }
    # print current plot if print_plot == T
    if (print_plot == T) print(p)
    # append all plot
    all_p[[my_title]] <- p 
  }
  
  # return vector of plot(s)
  if (length(all_p) == 1) all_p <- all_p[[1]]
  if (return_plots ) return(all_p)
}


# Analysis wrapper ----

# to do 

# Genetic analysis functions ----

require(here)

aa_convert <- function(x){
  amino_acid_code_df <- read_csv(paste0(here::here("Functions/amino-acid-code.csv")))
  aa_code_vec <- amino_acid_code_df$`One letter code`
  names( aa_code_vec ) <- str_to_title( amino_acid_code_df$`Three letter code`)
  
  for (i in 1:length( aa_code_vec ) ){
    x <- str_replace_all( x, names( aa_code_vec )[i], aa_code_vec[i] )
  }
  return(x)
}

# given a dataframe of isolates in iso1-iso2 format, generate a new dataframe with symmetric pairs and pair ids

create_gen_pairs <- function(df, # a dataframe with two columns: iso1 and iso2
                             initial_id_number = 1, # numeric with the first id number. This is usually 1, but will be higher if we are adding new pairs
                             write_file = NULL){# character with full path of ".tab" file to write output to. Useful to transfer the dataframe to the server
  
  df1 <- df
  
  df2 <- df %>%
    select(iso1 = iso2, iso2 = iso1)
  
  df_isolates <- bind_rows(df1, df2) %>%
    arrange(iso1) %>%
    mutate(PAIR_ID = str_c("GP-",
                           formatC(row_number() + initial_id_number - 1,
                                   digits = 3,
                                   format = "d",
                                   flag = "0"))) %>%
    relocate(PAIR_ID)
  
  if (!is.null(write_file)){
    
    df_isolates %>%
      write_tsv(write_file, col_names = F)
  }
  
  return(df_isolates)
}