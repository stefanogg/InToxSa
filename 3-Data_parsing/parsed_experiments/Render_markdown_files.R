# Run parsing markdown documents

# Library
library(rmarkdown)

# render all files
for (file in list.files("Data_parsing/parsed_experiments/", pattern = "GP", recursive = T, full.names = T) %>%
     str_subset(".Rmd")) render(file, output_format = "html_document")
