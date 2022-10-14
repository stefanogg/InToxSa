# Prevalent of ST in the genetic pairs

# Library
library(tidyverse)

# Data
data <- read_csv("Genetic_pairs_analysis_Oct_2020/metadata/genetic_pairs_id_metadata.csv") 

data <- data %>%
  select(iso1_ST, iso1) %>%
  distinct()

# Counts
data %>%
  count(iso1_ST, sort = T) %>%
  mutate(frac = scales::percent(n/sum(n)))
