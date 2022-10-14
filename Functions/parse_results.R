# A function to parse Excel tables

# Library
library(tidyverse)
library(magrittr)
library(readxl)

# Functions

# Parse PI kinetics ----

parse_kinetics <- function(file, # path to Excel file
                           n_skip = 14, # number of lines to skip (includes Excel column names)
                           interval = 6, # numeric with interval between measurements (in minutes)
                           duration, # numeric vector with duration of 
                           # the experiment in hours, minutes (e.g. 22 hours and 12 minutes: 22, 12)
                           wavelength = 535,# wavelengths used (numeric vector),
                           exp,
                           cell_number,
                           plate,
                           plot_data = TRUE
                           ){ 
  
  # Timepoint calculation in minutes
  duration_h <- duration[1]*60 # hour component of duration, transformed in minutes
  duration_min <- duration[2] # minute component of duration
  timepoint_max <- duration_h + duration_min
  # generates vector of timepoints (in minutes) of fluorescence measurement. Here: every 6 minutes until 22 hours and 18 minutes
  timepoints <- seq(0, timepoint_max, interval) 
  # generate vector of column names based on timepoints and wavelength. Here, both wave length were used
  col_names_WL <- sapply(wavelength, function(x){
    str_c("WL", x, "_", timepoints)
  })
  col_names <- c(
    "well",
    "sample",
    col_names_WL
  )
  data <- read_xlsx(path = file, col_names = col_names, skip = n_skip) 
  
  # transform into long format and add new variables
  data <- data %>%
      gather(key = condition,
             value = f_signal,
             starts_with("WL")) %>%
      separate(condition, into = c("wave_length", "time"), remove = FALSE) %>%
      filter(wave_length == "WL535") %>% # filter only used wavelength
      mutate(timepoint = as.numeric(time),
             well_row = str_extract(well, "[A-H]"),
             well_col = str_extract(well, "\\d{2}"),
             experiment = exp,
             cell_number = cell_number,
             plate = plate,
             well_id = paste(experiment, plate, well, sep = "_")) %>% # added a variable that identify unique experimental reading
      select(-c(condition, time, sample, wave_length))
  
  if (plot_data){
    p <- data %>%
      ggplot(aes(x = timepoint, y = f_signal)) +
      geom_point(size = .5, color = "blue") +
      facet_grid(well_row ~ well_col) +
      #scale_colour_manual(values = c("orange", "blue")) +
      theme_bw()
    print(p)
  }
  return(data)
}

# Parse LDH results ----

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

# Parse growth curves ----

parse_growth <- function(file,
                     n_skip = 14,
                     interval = 15,
                     duration,
                     normalised = FALSE,
                     exp,
                     plate){
  
  # fix arguments for parsing
  dir  <- getwd()
  path <- str_c(dir,
                file)
  
  duration_h <- duration[1]*60 # hour component of duration, transformed in minutes
  duration_min <- duration[2] # minute component of duration
  timepoint_max <- duration_h + duration_min
  # generates vector of timepoints (in minutes) of OD measurement.
  timepoints <- seq(0, timepoint_max, interval) 
  col_names_OD <- str_c("OD_", timepoints)
  col_names <- c("well_row", "well_col", "sample", col_names_OD)
  
  # parse data
  data <- read_xlsx(file, skip = n_skip, col_names = col_names) %>%
    pivot_longer(cols = starts_with("OD"), 
                 names_to = "time", 
                 values_to = "OD") %>%
    mutate(well = str_c(well_row, 
                        formatC(well_col, 
                                width = 2, 
                                format = "d", 
                                flag = "0")),
           time = str_remove(time, "OD_"),
           OD_blank_normalised = normalised,
           experiment = exp,
           plate = plate) %>%
    select(-c(well_row, well_col, sample)) 
  
  return(data)
  
}

# Parse inocula OD ----

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