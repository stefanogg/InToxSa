# Function to plot plate results
# df is a dataframe in long format with time-dependent results of the cell toxicity assay

plot_toxicity <- function(df, plate_number, wave_length = "WL535", save_plot = FALSE){
  title <- str_c("Plate ", plate_number)
  df <- df %>%
    filter(plate == plate_number) %>%
    mutate(experiment = str_c(experiment, "\n(", cell_number, " cells)"))
  p <- df %>%
    ggplot(aes(x = timepoint, y = f_signal)) +
    geom_point(size = .5, color = "red") +
    growthcurve::stat_growthcurve(model = "spline", color = "black") +
    facet_grid(experiment ~ sample_id, scales = "free_y") +
    labs(title = title,
         x = "Time (min)",
         y = "Fluorescence intensity") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  if (save_plot) {
    file <- str_c("figures/PI_kinetics_plate", plate_number, ".pdf" )
    ggsave(file, width = 20, height = 8)
  }
  print(p)
}

plot_toxicity_facet_wrap <- function(df, plate_number, wave_length = "WL535", n_row = 1, save_plot = FALSE){
  title <- str_c("Plate ", plate_number)
  df <- df %>%
    filter(plate == plate_number) %>%
    mutate(experiment = str_c(experiment, "\n(", cell_number, " cells)"))
  df <- df[df$wave_length == wave_length,]
  p <- df %>%
    ggplot(aes(x = timepoint, y = f_signal)) +
    geom_point(size = .5, color = "red") +
    growthcurve::stat_growthcurve(model = "spline", color = "black") +
    facet_wrap(strain_group ~ sample_id, nrow = n_row) +
    scale_x_continuous(breaks = c(0,10,20)) +
    labs(title = title,
         x = "Time (min)",
         y = "Fluorescence intensity") +
    theme_bw() 
  if (save_plot) {
    file <- str_c("figures/PI_kinetics_plate", plate_number, ".pdf" )
    ggsave(file, width = 20, height = 8)
  }
  print(p)
}
