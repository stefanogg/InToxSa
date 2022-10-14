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


