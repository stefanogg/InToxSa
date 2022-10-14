# Here we plot lineage effects for cytoxicity in VANANZ

library(tidyverse)
rm(list = ls())
setwd(str_c(here::here(), "/InToxSa_paper/GWAS/Lineages"))
getwd()

# Import raw data from pyseer
f <- "~/Documents/Transfer_with_server/lineage_effects.txt"
dir <- "raw_data/"
file.copy(f, dir)
df_lineage_eff <- read_tsv(str_c(dir, basename(f))) %>%
  rename(p_value = `p-value`)

# Plot
# get order
df_lineage_eff %>%
  ggplot(aes(x = fct_rev(fct_reorder(lineage, p_value)),
             y = -log10(p_value))) +
 geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(.05/10), colour = "red", linetype = "dotted") +
  labs(x = "Principal component", 
        y = "-log10(p)") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(text = element_text(face = "bold"))

dir <- "figures/"
dir.create(dir)
ggsave(str_c(dir, "lineage_effects_barplot.pdf"), width = 4.5, height = 2.9)

p <- df_lineage_eff %>%
  mutate(lineage = str_replace(lineage, "MDS", "PC")) %>%
  ggplot(aes(x = fct_reorder(lineage, p_value),
             y = -log10(p_value))) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(.05/10), colour = "red", linetype = "dotted") +
  labs(x = "Principal component", y = "-log10(p)") +
  theme_bw(base_size = 11) +
  theme(text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))
p

saveRDS(p, "processed_data/lineage_effects_barplot.Rda")

df_lineage_eff %>%
  arrange(p_value) %>%
  write_csv("processed_data/lineage_effects_data.csv")
