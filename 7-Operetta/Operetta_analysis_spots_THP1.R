
library(tidyverse)

setwd("/home/rguerillot/Documents/Travail/Abdou_project/Staph_infection_project/github_analysis/VANANZ_phenotypes")

# import cell count operetta data
read_tsv_filename <- function(flnm) {
  read_tsv(flnm, skip = 9, trim_ws = T) %>% 
    mutate(filename = flnm)
}

Spot_count.df <-list.files(path = "Operetta/raw_data/210921 THP1 n2/", pattern = "*Spots", full.names = T, recursive = T) %>% 
  map_df(~read_tsv_filename(.)) 

unique(Spot_count.df$filename)

# Tidy Operetta combined dataframe

Spot_count_clean.df <- Spot_count.df %>%
  separate(col = filename, sep = "/", into = c("a","b","c", "d", "experiment_id")) %>%
  separate(col = experiment_id, sep = "_", into = c("plate_date", "cell_type", "timepoint", "e", "f"), remove = F) %>%
  select(-a, -b, -d, -e, -f, -Timepoint) %>%
  mutate(Row = as.character(Row)) %>%
  mutate(Strain = ifelse(Strain == "Non-infected", yes = "non-infected", no = as.character(Strain))) %>%
  mutate(`Cell Type` = ifelse(`Cell Type` == "THP1 Casp1-/-", yes = "THP1 casp1-/-", no = as.character(`Cell Type`))) 
  

raw_to_ABC.df <- data_frame(Row = c(1,2,3,4,5,6,7,8), row = c("A","B","C","D","E","F","G","H")) %>%
  mutate(Row = as.character(Row))

Spot_count_clean.df <- left_join(Spot_count_clean.df, raw_to_ABC.df) %>%
  select(-Row) %>%
  select(Row = row, Column, everything()) %>%
  mutate(Column = str_pad(Column, 2, pad = "0")) %>%
  mutate(Well = paste0(Row, Column)) %>%
  select(Well, everything(), -Replicate) %>%
  mutate(sample_id = paste(experiment_id, Well, `Cell Type`, Strain, sep = "#"))
  
# create a df of replicates and merge with clean data 
sample_replicate_df <- Spot_count_clean.df %>%
  select(experiment_id, Well, `Cell Type`, Strain) %>%
  distinct() %>%
  group_by(experiment_id, `Cell Type`, Strain) %>%
  mutate(replicate = row_number()) %>%
  ungroup()

Spot_count_clean.df <- merge(Spot_count_clean.df, sample_replicate_df, by = c("experiment_id", "Well", "Cell Type", "Strain")) %>%
  mutate(timepoint = factor(timepoint, levels = c("1.5h", "5h", "24h"))) %>%
  mutate(Strain = factor(Strain, levels = c("WT", "agrA", "non-infected"))) %>%
  mutate(`Cell Type` = factor(`Cell Type`, levels = c("THP1 (Cas9)", "THP1 casp1-/-", "THP1 casp-4/5 -/-"))) 

# count bacteria per field
bac_count_field <- Spot_count_clean.df  %>% group_by(timepoint,  Well, `Cell Type`, Strain, replicate, Field) %>%
  count(name = "number of bacteria/field")

ggplot(bac_count_field, aes(x = Strain, y = `number of bacteria/field`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  facet_grid(~timepoint)

# count infected cells
bac_count_infected <- Spot_count_clean.df  %>% group_by(timepoint,  Well, `Cell Type`, Strain, replicate, Field) %>%
  distinct(`Spots - Object No in Cells`) %>%
  count(name = "number of infected cells/field")

# count bacteria per cell
bac_count_cell <- Spot_count_clean.df  %>% 
  group_by(timepoint,  Well, `Cell Type`, Strain, replicate, Field, `Spots - Object No in Cells`) %>%
  count(name = "number of bacteria/infected cell") %>%
  ungroup() %>%
  group_by(timepoint,  Well, `Cell Type`, Strain, replicate, Field, `number of bacteria/infected cell`) %>%
  count(name = "number of infected cells")
  

ggplot(bac_count_cell, aes(x = `number of bacteria/infected cell`, y = `number of infected cells`, fill = `Cell Type`,
                           group =  interaction(`Cell Type`, `number of bacteria/infected cell`))) +
  geom_boxplot(outlier.shape = NA, position = "dodge2")+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = .15),alpha=0.3, , colour = "black") +
  facet_grid(Strain ~timepoint)+
  xlim(0,7)

ggplot(bac_count_cell, aes(x = `number of bacteria/infected cell`, y = `number of infected cells`, fill = `Strain`,
                           group =  interaction(`Strain`, `number of bacteria/infected cell`))) +
  geom_boxplot(outlier.shape = NA, position = "dodge2")+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = .15),alpha=0.3, , colour = "black") +
  facet_grid(`Cell Type` ~timepoint)+
  xlim(0,7)

# Plot as specified by Abdou

p1 <- ggplot(bac_count_field, aes(x = timepoint, y = `number of bacteria/field`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  ylab("intracellular bacteria per field") +
  theme_bw() 
p1
ggsave("p1.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p1.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)

p2 <- ggplot(bac_count_field %>%
               filter(`Cell Type` != "THP1 casp-4/5 -/-")
             , aes(x = timepoint, y = `number of bacteria/field`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  ylab("intracellular bacteria per field") +
  theme_bw()
p2
ggsave("p2.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p2.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)

p3 <- ggplot(bac_count_field %>%
               filter(`Cell Type` != "THP1 casp1-/-")
             , aes(x = timepoint, y = `number of bacteria/field`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  ylab("intracellular bacteria per field") +
  theme_bw()
p3
ggsave("p3.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p3.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)


p4 <- ggplot(bac_count_cell, aes(x = `number of bacteria/infected cell`, y = `number of infected cells`, fill = `Cell Type`,
                           group =  interaction(`Cell Type`, `number of bacteria/infected cell`))) +
  geom_boxplot(outlier.shape = NA, position = "dodge2")+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = .15),alpha=0.3, , colour = "black") +
  facet_grid(Strain~timepoint) +
  theme_bw() +
  scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10))
p4
ggsave("p4.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 14, height = 6)
ggsave("p4.png", path = "Operetta/Operetta_THP1_all_plots/", width = 14, height = 6)


p5 <- ggplot(bac_count_cell %>%
               filter(`Cell Type` != "THP1 casp-4/5 -/-")
             , aes(x = `number of bacteria/infected cell`, y = `number of infected cells`, fill = `Cell Type`,
                                 group =  interaction(`Cell Type`, `number of bacteria/infected cell`))) +
  geom_boxplot(outlier.shape = NA, position = "dodge2")+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = .15),alpha=0.3, , colour = "black") +
  facet_grid(Strain~timepoint) +
  theme_bw() +
  scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10))
p5
ggsave("p5.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 14, height = 6)
ggsave("p5.png", path = "Operetta/Operetta_THP1_all_plots/", width = 14, height = 6)


p6 <- ggplot(bac_count_cell %>%
               filter(`Cell Type` != "THP1 casp1-/-")
             , aes(x = `number of bacteria/infected cell`, y = `number of infected cells`, fill = `Cell Type`,
                   group =  interaction(`Cell Type`, `number of bacteria/infected cell`))) +
  geom_boxplot(outlier.shape = NA, position = "dodge2")+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = .15),alpha=0.3, , colour = "black") +
  facet_grid(Strain~timepoint) +
  theme_bw() +
  scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10))
p6
ggsave("p6.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 14, height = 6)
ggsave("p6.png", path = "Operetta/Operetta_THP1_all_plots/", width = 14, height = 6)

# idem without points

p7 <- ggplot(bac_count_cell, aes(x = `number of bacteria/infected cell`, y = `number of infected cells`, fill = `Cell Type`,
                                 group =  interaction(`Cell Type`, `number of bacteria/infected cell`))) +
  geom_boxplot(outlier.shape = NA, position = "dodge2")+
  #geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = .15),alpha=0.3, , colour = "black") +
  facet_grid(Strain~timepoint) +
  theme_bw() +
  scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10))
p7
ggsave("p7.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 14, height = 6)
ggsave("p7.png", path = "Operetta/Operetta_THP1_all_plots/", width = 14, height = 6)

p8 <- ggplot(bac_count_cell %>%
               filter(`Cell Type` != "THP1 casp-4/5 -/-")
             , aes(x = `number of bacteria/infected cell`, y = `number of infected cells`, fill = `Cell Type`,
                   group =  interaction(`Cell Type`, `number of bacteria/infected cell`))) +
  geom_boxplot(outlier.shape = NA, position = "dodge2")+
  #geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = .15),alpha=0.3, , colour = "black") +
  facet_grid(Strain~timepoint) +
  theme_bw() +
  scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10))
p8
ggsave("p8.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 14, height = 6)
ggsave("p8.png", path = "Operetta/Operetta_THP1_all_plots/", width = 14, height = 6)

p9 <- ggplot(bac_count_cell %>%
               filter(`Cell Type` != "THP1 casp1-/-")
             , aes(x = `number of bacteria/infected cell`, y = `number of infected cells`, fill = `Cell Type`,
                   group =  interaction(`Cell Type`, `number of bacteria/infected cell`))) +
  geom_boxplot(outlier.shape = NA, position = "dodge2")+
  #geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = .15),alpha=0.3, , colour = "black") +
  facet_grid(Strain~timepoint) +
  theme_bw() +
  scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10))
p9
ggsave("p9.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 14, height = 6)
ggsave("p9.png", path = "Operetta/Operetta_THP1_all_plots/", width = 14, height = 6)

# Plot total number of infected cells

p10 <- ggplot(bac_count_infected, aes(x = timepoint, y = `number of infected cells/field`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  ylab("infected cells per field") +
  theme_bw()
p10
ggsave("p10.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p10.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)

p11 <- ggplot(bac_count_infected %>% 
                filter(`Cell Type` != "THP1 casp-4/5 -/-"), aes(x = timepoint, y = `number of infected cells/field`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  ylab("infected cells per field") +
  theme_bw()
p11
ggsave("p11.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p11.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)


p12 <- ggplot(bac_count_infected %>% 
                filter(`Cell Type` != "THP1 casp1-/-"), aes(x = timepoint, y = `number of infected cells/field`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  ylab("infected cells per field") +
  theme_bw()
p12
ggsave("p12.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p12.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)


# note: Number of infected increase for agrA between 5h and 24h -> less seeded cells on 5h plate?? => need to correct by total number of cells

# Plot total nb cells and pct of infected cells ----

Cell_count.df <-list.files(path = "Operetta/raw_data/210921 THP1 n2/", pattern = "*Cells", full.names = T, recursive = T) %>% 
  map_df(~read_tsv_filename(.)) 

unique(Spot_count.df$filename)

Cell_count_clean.df <- Cell_count.df %>%
  separate(col = filename, sep = "/", into = c("a","b","c", "d", "experiment_id")) %>%
  separate(col = experiment_id, sep = "_", into = c("plate_date", "cell_type", "timepoint", "e", "f"), remove = F) %>%
  select(-a, -b, -c, -d, -e, -f, -Timepoint) %>%
  mutate(Row = as.character(Row)) %>%
  mutate(Strain = ifelse(Strain == "Non-infected", yes = "non-infected", no = as.character(Strain))) %>%
  mutate(`Cell Type` = ifelse(`Cell Type` == "THP1 Casp1-/-", yes = "THP1 casp1-/-", no = as.character(`Cell Type`))) %>%
  mutate( Cells = ifelse(`Cells - Number of Spots` > 0, yes = "infected cells", "non-infected cells"))  %>%
  mutate(timepoint = factor(timepoint, levels = c("1.5h", "5h", "24h"))) %>%
  mutate(Strain = factor(Strain, levels = c("WT", "agrA", "non-infected"))) %>%
  mutate(`Cell Type` = factor(`Cell Type`, levels = c("THP1 (Cas9)", "THP1 casp1-/-", "THP1 casp-4/5 -/-"))) 

total_cells_per_field <- Cell_count_clean.df  %>% group_by(timepoint,  Row, Column, `Cell Type`, Strain, Replicate, Field) %>%
  count(name = "number of cells") %>%
  ungroup()

infected_cells_per_field <- Cell_count_clean.df %>%
  filter(Cells == "infected cells") %>%
  group_by(timepoint,  Row, Column, `Cell Type`, Strain, Replicate, Field) %>%
  count(name = "number of infected cells") %>%
  ungroup()

pct_infected_cells_per_field <- merge(total_cells_per_field, infected_cells_per_field,
                                      by = c("timepoint", "Row", "Column", "Cell Type", "Strain", "Replicate", "Field"),
                                      all.x = T) %>%
  mutate(`number of infected cells` = ifelse(is.na(`number of infected cells`),yes = 0, no = `number of infected cells`)) %>%
  mutate(`% of infected cells` = (`number of infected cells`/`number of cells`)*100)


p13 <- ggplot(pct_infected_cells_per_field, aes(x = timepoint, y = `number of cells`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  #  ylab("intracellular bacteria per field") +
  theme_bw()
p13
ggsave("p13.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p13.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)


p14 <- ggplot(pct_infected_cells_per_field, aes(x = timepoint, y = `number of infected cells`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  #  ylab("intracellular bacteria per field") +
  theme_bw()
p14
ggsave("p14.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p14.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)


p15 <- ggplot(pct_infected_cells_per_field, aes(x = timepoint, y = `% of infected cells`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
#  ylab("intracellular bacteria per field") +
  theme_bw()
p15
ggsave("p15.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p15.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)

p16 <- ggplot(pct_infected_cells_per_field %>%
                filter(`Cell Type` != "THP1 casp-4/5 -/-"), aes(x = timepoint, y = `number of cells`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  #  ylab("intracellular bacteria per field") +
  theme_bw()
p16
ggsave("p16.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p16.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)

p17 <- ggplot(pct_infected_cells_per_field %>%
                filter(`Cell Type` != "THP1 casp-4/5 -/-"), aes(x = timepoint, y = `number of infected cells`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  #  ylab("intracellular bacteria per field") +
  theme_bw()
p17
ggsave("p17.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p17.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)

p18 <- ggplot(pct_infected_cells_per_field %>%
                filter(`Cell Type` != "THP1 casp-4/5 -/-"), aes(x = timepoint, y = `% of infected cells`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  #  ylab("intracellular bacteria per field") +
  theme_bw()
p18
ggsave("p18.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p18.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)

p19 <- ggplot(pct_infected_cells_per_field %>%
                filter(`Cell Type` != "THP1 casp1-/-"), aes(x = timepoint, y = `number of cells`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  #  ylab("intracellular bacteria per field") +
  theme_bw()
p19
ggsave("p19.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p19.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)

p20 <- ggplot(pct_infected_cells_per_field %>%
                filter(`Cell Type` != "THP1 casp1-/-"), aes(x = timepoint, y = `number of infected cells`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  #  ylab("intracellular bacteria per field") +
  theme_bw()
p20
ggsave("p20.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p20.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)

p21 <- ggplot(pct_infected_cells_per_field %>%
                filter(`Cell Type` != "THP1 casp1-/-"), aes(x = timepoint, y = `% of infected cells`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~Strain) +
  #  ylab("intracellular bacteria per field") +
  theme_bw()
p21
ggsave("p21.pdf", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)
ggsave("p21.png", path = "Operetta/Operetta_THP1_all_plots/", width = 10, height = 4)




# Check signals differences infected vs non-infected and try to remove background bact. counts in non-infected ----
ninf.df <- Spot_count_clean.df %>%
  filter(Strain == "non-infected") %>%
  select(starts_with("Spots")) %>%
  gather() %>%
  mutate(Strain = "non-infected")

inf.df <- Spot_count_clean.df %>%
  filter(Strain != "non-infected") %>%
  select(starts_with("Spots")) %>%
  gather() %>%
  mutate(Strain = "infected")

ninf_inf.df <- rbind(ninf.df, inf.df) %>%
  filter(!is.na(value))

unique(ninf_inf.df$key)

ggplot(ninf_inf.df, aes(x = value, colour = Strain)) +
  geom_density() +
  facet_wrap(~ key, scales = "free") 

ggplot(ninf_inf.df, aes(x = value, colour = Strain)) +
  geom_density() +
  facet_wrap(~ key, scales = "free") +
  xlim(0, 100)

# based on distribution background spot/bacteria could be reduced by removing Spot area > 37.5

Spot_count_clean_clean.df <-  Spot_count_clean.df %>%
  filter(`Spots - Spot Area [px²]` < 37.5)

# Count and re plot after filter
bac_count_field <- Spot_count_clean_clean.df  %>% group_by(timepoint,  Well, `Cell Type`, Strain, replicate, Field) %>%
  count(name = "number of bacteria/field")

ggplot(bac_count_field, aes(x = Strain, y = `number of bacteria/field`, fill = `Cell Type`)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = 0),alpha=0.3) +
  #  geom_jitter(width = .2, alpha= .3) +
  facet_grid(~timepoint)

# count bacteria per cell
bac_count_cell <- Spot_count_clean_clean.df  %>% 
  group_by(timepoint,  Well, `Cell Type`, Strain, replicate, Field, `Spots - Object No in Cells`) %>%
  count(name = "number of bacteria/infected cell") %>%
  ungroup() %>%
  group_by(timepoint,  Well, `Cell Type`, Strain, replicate, Field, `number of bacteria/infected cell`) %>%
  count()


ggplot(bac_count_cell, aes(x = `number of bacteria/infected cell`, y = n, fill = `Cell Type`,
                           group =  interaction(`Cell Type`, `number of bacteria/infected cell`))) +
  geom_boxplot(outlier.shape = NA, position = "dodge2")+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = .15),alpha=0.3, , colour = "black") +
  facet_grid(Strain ~timepoint)+
  xlim(0,7)

ggplot(bac_count_cell, aes(x = `number of bacteria/infected cell`, y = n, fill = `Strain`,
                           group =  interaction(`Strain`, `number of bacteria/infected cell`))) +
  geom_boxplot(outlier.shape = NA, position = "dodge2")+
  geom_point(position=position_jitterdodge(jitter.width = .15, jitter.height = .15),alpha=0.3, , colour = "black") +
  facet_grid(`Cell Type` ~timepoint)


# Note: doesn't work remove most of the data -> need to rework background removal during image processing steps

