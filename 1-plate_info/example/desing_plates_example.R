# Here we provide an example of plate design and randomisation

# We start with a dataframe of isolates (that may be grouped), a list of controls. We assume 3 replicates for the strains, but the number of replicates 
# of controls can vary

# Source functions
source("plates_design_functions.R")

# Raw data and variables
df_isolates <- read_csv("df_isolates.csv")
controls <- c("JE2", "TOX-4", "Complete lysis", "Non infected")
n_replicates_controls <- c(3, 3, 3, 3) # a numeric of the same length as the control vector for the number of replicates for each control

# Distribute strains across plates. This function returns a dataframe of the same length of df_isolates
# with a supplementary column indicating the plate assignment (as plate number). If plot = T (the default) it prints a plot to check the distribution
df_plate_distribution <- plate_distributor(df_isolates = df_isolates,
                  n_replicates_controls = n_replicates_controls,
                  plate_number_format = "PLATE-")

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
    fixed = NULL,
    controls = controls,
    n_replicates_controls = n_replicates_controls
  ) %>%
    mutate(plate_number = x)
  
  write_csv(well_info,
            str_c(dir, "well_info_plate", x, ".csv"))
  
  return(well_info)
}) %>%
  bind_rows()

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

