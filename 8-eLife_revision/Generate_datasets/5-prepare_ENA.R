# Prepare input files for ENA submission

# Library 
library(tidyverse)

setwd(stringr::str_c("~/Documents/Github/InToxSa", "/6-InToxSa_paper/Generate_datasets"))

rm(list = ls())

# Sequences dataset (allelic exchange mutants)
df_sequences <- readxl::read_xlsx("processed_data/supp_table_4/supp_table_4.xlsx", sheet = "allelic exchange mutants") %>%
  filter(str_detect(sample_id, "BPH3370-")) %>%
  select(sample_id, strain_description = mutation)

# Construct ENA dataset
df_ena_sample <- df_sequences %>%
  transmute(tax_id = 1280,
            scientific_name = "Staphylococcus aureus",
            sample_alias = sample_id,
            sample_title = "InToxSa allelic exchange",
            sample_description = strain_description,
            isolation_source = "Laboratory",
            collection_date = 2022,
            `geographic location (country and/or sea)` = "Australia",
            `host health state` = "diseased",
            `host scientific name` = "Homo sapiens",
            isolate = sample_id)

# save 
# df_ena_sample %>%
#   write_tsv("~/OneDrive - The University of Melbourne/Postdoc_work/VANANZ_phenotypes/ENA_submission/Samples_submission/Checklist_ENA-prokaryotic pathogen minimal sample checklist_1670906494278.tsv", append = T, col_names = F)

# Construct reads dataset
df_ena_accession <- read_csv("~/OneDrive - The University of Melbourne/Postdoc_work/VANANZ_phenotypes/ENA_submission/Samples_submission/Reports/samples-2022-12-20T00_24_22.csv") 
df_ena_reads <- df_ena_accession %>%
  transmute(alias,
            sample = id,
            study = "PRJEB58038",
            instrument_model = "NextSeq 500",
            library_name = "unspecified",
            library_soure = "GENOMIC",
            library_selection = "RANDOM",
            library_strategy = "WGS",
            library_layout = "PAIRED")

# Import files names and md5
f <- "~/Documents/Transfer_with_server/isolates.md5.tab"
dir <- "raw_data/"
file.copy(f, dir)
df_ena_files <- read_tsv(str_c(dir, basename(f)),
                         col_names = c("alias", "R1", "R2")) %>%
  separate(R1, into = c("forward_file_md5", "forward_file_name"), sep = "  ") %>%
  separate(R2, into = c("reverse_file_md5", "reverse_file_name"), sep = "  ") %>%
  mutate(across(.cols = ends_with("name"),
                .fns = ~basename(.))) %>%
  select(alias, forward_file_name, forward_file_md5, reverse_file_name, reverse_file_md5)

# Merge
df_ena_reads <- df_ena_reads %>%
  inner_join(df_ena_files) %>%
  select(-alias) %>%
  write_tsv("~/OneDrive - The University of Melbourne/Postdoc_work/VANANZ_phenotypes/ENA_submission/Reads_submission/fastq2_template_1671495526442.tsv", append = T, col_names = F)

# Check accession numbers
f <- "~/OneDrive - The University of Melbourne/Postdoc_work/VANANZ_phenotypes/ENA_submission/Reads_submission/Reports/runs-2022-12-20T00_34_39.csv"
file.copy(f, dir)
df_ena_run <- read_csv(str_c(dir, basename(f))) 

# Construct a database with reads metadata and accession numbers
df_ena <- df_ena_sample %>%
  inner_join(df_ena_accession, by = c("sample_alias" = "alias")) %>%
  rename(sample_accession2 = id, sample_accession = secondaryId) %>%
  inner_join(df_ena_reads, by = c("sample_accession2" = "sample")) %>%
  inner_join(df_ena_run %>% select(sampleId, run_accession = id, alias, firstCreated_run = firstCreated, firstPublic_run = firstPublic, studyId, experimentId), by = c("sample_accession2" = "sampleId"))


# save
df_ena %>%
  write_csv("processed_data/ena_accession_metadata.csv")
