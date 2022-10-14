get_parameters_stat <- function(df_param,
                                ref = "JE2",
                                test = "wilcox.test",
                                group = NULL) { # a character vector containing the name of grouping variables (eg. c("experiment"))
  require(ggpubr)
  t <- rbind(compare_means(formula =  max_rate_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  time_of_max_rate_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  min_rate_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  time_of_min_rate_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  max_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  time_of_max_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  min_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  time_of_min_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group),
             compare_means(formula =  AUC_death ~ sample_id, 
                           data = df_param, 
                           ref.group = ref, 
                           method = test,
                           group.by = group))
  colnames(t)[1] <- "tested_parameter"
  colnames(t)[2] <- "reference"
  colnames(t)[3] <- "sample_id"
  
  if (!is.null(group)){
    colnames(t)[1] <- "grouping"
    colnames(t)[2] <- "tested_parameter"
    colnames(t)[3] <- "reference"
    colnames(t)[4] <- "sample_id"
    
  }
  return(t)
}

t <- parameters_stat_NEBRASKA.df <- parameters_PI.df %>%
  filter(strain_group == "NEBRASKA" | sample_id == "JE2") %>%
  compare_means(formula =  AUC_death ~ sample_id, 
                data = ., 
                ref.group = "JE2", 
                method = "wilcox.test",
                group.by = NULL)
t_signif <- t %>% filter(p.format < 0.01) %>% .$group2 
#write.csv(t_signif, "Documents/Travail/plot_grant_2020_NEB.csv", row.names = F)

ggboxplot(data = parameters_PI.df %>% filter(sample_id %in% c(t_signif,"JE2"))%>%
          filter(sample_id != "Maximum cell death"), 
          y = "AUC_death", 
          x = "gene_name", 
          color = "strain_group") +
          #geom_point(aes(colour = experiment)) +
  stat_compare_means(ref.group = "JE2", 
                     method = "wilcox.test",
                     label = "p.signif")

a <- filter(strain_info_df, sample_id %in% t_signif) %>%
  select(sample_id,gene_name)


d2 <-read.csv( "Documents/Travail/plot_grant_2020_NEB.csv")
dat <- parameters_PI.df %>% filter(sample_id %in% c(t_signif,"JE2"))%>%
  merge(., d2, by.x= "sample_id", by.y = "group2",all.x = T, all.y = T)


ggboxplot(data = dat, order = c("agrC","spa","dsbA","fmtC","hp1","hp2","hp3","WT"), 
          legend = "none",
          y = "AUC_death", 
          x = "gene",
          outlier.shape=NA,
          color = "strain_group") +
  #geom_point(aes(colour = experiment)) +
  stat_compare_means(ref.group = "WT", 
                     method = "wilcox.test",
                     label = "p.signif")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),legend.position = "none" )+
  ylab("cytotox.") +
  xlab("")
  
  

