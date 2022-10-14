# Generate dataframes of metadata (annotation) for pyseer output

# Mutations: metadata from the vcf file with all non-syn mutations + FPR3757 homolog

# Genes: summarise vcfs for the different genotypes + FPR3757 homolog

# Library
library(tidyverse)
library(vcfR)
rm(list = ls())

source("~/OneDrive - The University of Melbourne//R/Functions/aa_convert.R")

setwd(str_c(here::here(), "/GWAS/pyseer_output"))
dir_raw <- "raw_data/"
dir.create(dir_raw)
dir <- str_c(dir_raw, "annotation/")
dir.create(dir)

# Copy vcf file
f <- "~/Documents/Transfer_with_server/nosyn.core90.vcf.gz"
file.copy(f, dir, overwrite = F)


# Parse vcf
vcf <- read.vcfR(str_c(dir, basename(f)))


df_vcf <- vcfR2tidy(vcf, info_only = T)
df_mutations <- df_vcf$fix 
df_meta <- df_vcf$meta
rm(df_vcf)

df_mutations_metadata <- df_mutations %>%
  distinct(CHROM, POS, REF, ALT, TYPE, ANN, AF) %>%
  mutate(EFFTYPE = str_split_fixed(ANN, pattern = "\\|", n = 16)[,2],
         LOCUS_TAG = str_split_fixed(ANN, pattern = "\\|", n = 16)[,4],
         MUTATION = str_split_fixed(ANN, pattern = "\\|", n = 16)[,11]) %>%
  mutate(MUTATION_SHORT = aa_convert(str_remove(MUTATION, "p."))) %>%
  mutate(mutation_id = str_c(CHROM, POS, REF, ALT, sep = "_")) %>%
  mutate(AF = as.numeric(AF))

# add genes genotype tags using the same criteria as in bcftools
df_mutations_metadata <- df_mutations_metadata %>%
  mutate(rare_mutation = AF <= 0.01,
         truncation = str_detect(EFFTYPE, "frameshift_variant") |
                          str_detect(EFFTYPE, "stop_gained") |
                          str_detect(EFFTYPE, "stop_lost") |
                          str_detect(EFFTYPE, "start_lost"))

# df_mutations_wide <- df_mutations_metadata %>%
#   left_join(df_mutations) %>%
#   select(CHROM, POS, REF, ALT, TYPE, MUTATION_SHORT, LOCUS_TAG, AF, sample_id, gt_GT) %>%
#   pivot_wider(names_from = sample_id, values_from = gt_GT)

# Add homoplasy data
f <- "~/Documents/Transfer_with_server/consistencyIndexReport_24-05-22.txt"
file.copy(f, dir, overwrite = F)
df_homoplasy <- read_tsv(str_c(dir, basename(f))) %>%
  select(POS = End, consistency_index = ConsistencyIndex,
         mutation_counts = Counts, n_acquisitions = MinimumNumberChangesOnTree)

# Add neb annotation

df_blast_neb <- readRDS("~/OneDrive - The University of Melbourne/R/SAUREUS-GENERAL/ref_genomes/processed_data/BPH2947/df_blastp_BPH2947_neb_curated.Rda")



df_annotation <- df_mutations_metadata %>%
  left_join(df_homoplasy) %>%
  mutate(homoplasy = !is.na(consistency_index)) %>%
  relocate(homoplasy, .after = truncation) %>%
  left_join(df_blast_neb) %>%
  relocate(unique_gene_symbol)

# Save processed data
dir <- "processed_data/"
dir.create(dir)
dir <- str_c(dir, "annotation/")
dir.create(dir)

# df_annotation %>%
#   saveRDS(str_c(dir, "df_annotation.Rda"))

rm(df_meta, df_mutations, df_mutations_metadata, df_blast_neb, df_homoplasy)

# Extract genotype matrix rows = variants, cols = samples
geno <- extract.gt(vcf, return.alleles = F)
rownames(geno) <- df_annotation$mutation_id

# convert geno to rows = samples, cols = variants
str(geno)
geno2 <- t(geno)
str(geno2)
geno2[1:10, 1:10]

# Save processed data
geno2 %>%
  saveRDS(str_c(dir, "geno.Rda"))

rm(list = ls())

setwd(here::here())
