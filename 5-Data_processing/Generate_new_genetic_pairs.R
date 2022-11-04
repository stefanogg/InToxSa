# Short script to generate new genetic pairs on top of the 2400 that were generated initially

# Library
library(tidyverse)
source("Functions/all_functions.R")

# Useful code to print pairs for tribble()

# for (i in setdiff(missing_pairs, str_c(df_denovo_mutcounts$iso1, df_denovo_mutcounts$iso2, sep = "-"))){
#   string <- str_replace(i, '-', '\", \"' )
#   string <- str_c('\"', string, '\",\n')
#   cat(string)
# }

# PI pairs
# needs to add 
df_isolates <- tribble(
  ~iso1, ~iso2,
  "BPH2858", "BPH2720",
  "BPH2761", "BPH2752",
  "BPH2860", "BPH2784"
)
df1 <- df_isolates

df2 <- df_isolates %>%
  select(iso1 = iso2, iso2 = iso1)

df_isolates <- bind_rows(df1, df2) %>%
  arrange(iso1) %>%
  mutate(PAIR_ID = str_c("GP-",
                         formatC(row_number() + 3016,
                                 digits = 3,
                                 format = "d",
                                 flag = "0"))) %>%
  relocate(PAIR_ID)

# write_tsv(df_isolates,
#           "~/Documents/Transfer_with_server/genetic_pairs_GP-3017.tab",
#           col_names = F)

# One additional mortality pair (24 Feb 2021)
# We have noticed that BPH3678 (within the 11_BPH3676 cluster) was not included
df_isolates <- tribble(
  ~iso1, ~iso2,
  "BPH3678", "BPH3676"
)
df_isolates <- create_gen_pairs(df_isolates,
                 initial_id_number = 3023, write_file = "~/Documents/Transfer_with_server/genetic_pairs_GP-3023.tab")

# Add additional GC pairs
df_isolates <- tribble(
  ~iso1, ~iso2,
  "BPH2702", "BPH2816",
  "BPH2763", "BPH2823",
  "BPH2700", "BPH2876",
  "BPH3708", "BPH3541",
)
df_isolates <- create_gen_pairs(df_isolates,
                                initial_id_number = 3025, write_file = "~/Documents/Transfer_with_server/genetic_pairs_GP-3025.tab")

# Add additional mortality pairs
df_isolates <- tribble(
  ~iso1, ~iso2,
  "BPH2948", "BPH3455",
  "BPH2947", "BPH3455",
  "BPH3737", "BPH3328",
  "BPH3356", "BPH3328",
)
df_isolates <- create_gen_pairs(df_isolates,
                                initial_id_number = 3033, write_file = "~/Documents/Transfer_with_server/genetic_pairs_GP-3033.tab")
