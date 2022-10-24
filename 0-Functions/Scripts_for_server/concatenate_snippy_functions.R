# Function to remove duplicate mutations from concatenated snippy output
remove_duplicate_mutations <- function(file,
                                       pair_name = NULL,
                                       n_isolates = 2,
                                       write_file = F,
                                       compare_to_index = NULL) {
  require(tidyverse)
  require(magrittr)
  
  df <- read_tsv(file)
  
  # need a concatenated snippy output
  if (!is.null(compare_to_index)) {
    df_filtered <- df %>%
      group_by(CHROM, POS) %>%
      filter(!any(ISOLATE == index))
  }
  
  df_filtered <- df %>%
    group_by(CHROM, POS) %>%
    filter(n() < n_isolates) %>%
    # add pairname
    mutate(PAIR_ID = pair_name)
  
  
  
  
  if (write_file) {
    write_tsv(df_filtered,
              "all_snps_filtered.tab")
    return(NULL)
  }
  
  return(df_filtered)
  
}