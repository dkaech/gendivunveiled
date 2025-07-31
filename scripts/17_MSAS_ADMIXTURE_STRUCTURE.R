# Purpose: Analyse ADMIXTURE and STRUCTURE output of MSAS
# Author: Damian Käch <dkaech@ethz.ch>


# ADMIXTURE MSAS ----------------------------------------------------------

# Sample names file
sample_names_admix_ampse <- samples_names_ampse %>%
  filter(sample != "F-LT-1")

# List Q files
listed_files_focus_ampse <- list.files(path = paste0(path_admix_ampse, "focus_groups/"),
                                       pattern = ".Q$",
                                       full.names = T)

# Read listed files and set names
files_focus_ampse <- map(listed_files_focus_ampse, read.table)
names(files_focus_ampse) <- listed_files_focus_ampse %>%
  gsub(".*topweide_", "", .) %>%
  gsub("\\.Q", "", .) %>%
  gsub("\\.", "_", .)

# Adapt dfs in list
files_pure_ampse <- map(files_focus_ampse[c(na.omit(str_extract(names(files_focus_ampse),"pgh.*")))], 
                        function(df) df %>% 
                          mutate(names = str_extract(sample_names_admix_ampse$names, "^...-.*") %>% 
                                   na.omit() %>% 
                                   c()))

files_m1_ampse <- map(files_focus_ampse[c(na.omit(str_extract(names(files_focus_ampse),"m1_.*")))],
                      function(df) df %>%
                        mutate(names = str_extract(sample_names_admix_ampse$names, ".*(Ari|Rep).*") %>%
                                 na.omit() %>%
                                 c()))
files_m2_ampse <- map(files_focus_ampse[c(na.omit(str_extract(names(files_focus_ampse),"m2.*")))],
                      function(df) df %>%
                        mutate(names = str_extract(sample_names_admix_ampse$names, ".*(Ari|Art).*") %>%
                                 na.omit() %>%
                                 c()))
files_m3_ampse <- map(files_focus_ampse[c(na.omit(str_extract(names(files_focus_ampse),"m3.*")))],
                      function(df) df %>%
                        mutate(names = str_extract(sample_names_admix_ampse$names, ".*(Ari|Alg).*") %>%
                                 na.omit() %>%
                                 c()))

files_kcol_ampse <- c(files_pure_ampse, files_m1_ampse, files_m2_ampse, files_m3_ampse)
files_kcol_ampse <- map2(files_kcol_ampse, names(files_kcol_ampse),
                         function(df, df_name) df %>%
                           mutate(n_clusters = as.numeric(str_extract(df_name, "(?<=_).*")),
                                  # add column with information about the samples included in the ADMIXTURE analysis
                                  focus_group = str_extract(df_name, ".*(?=_)")) %>%
                           left_join(sample_names_admix, by = "names") %>%
                           pivot_longer(cols = contains("V"), names_to = "Cluster_name", values_to = "Cluster_value") %>%
                           mutate(Cluster_name = gsub("V", "Cluster_", Cluster_name) %>%
                                    factor(levels = paste0("Cluster_", 1:max(n_clusters)))))

# Bind all dfs
admix_ampse <- bind_rows(files_kcol_ampse) %>%
  mutate(names = case_when(str_detect(names, "LT-Rep50-3") ~ "LT-Rep50-1",
                           str_detect(names, "Rep-3") ~ "Rep-2",
                           str_detect(names, "Art25-3") ~ "Art25-1",
                           TRUE ~ names))

# Make sure to load admix_fct() and order of names (e.g., order_m1) from the GBS analysis

# M1
(admix_m1_ampse <- admix_ampse %>%
    filter(focus_group == "m1") %>%
    mutate(names = factor(names, level = order_m1)) %>%
    admix_fct(hide_legend = T))

# ggsave(plot = admix_m1_ampse, filename = paste0(path_figures, "MSAS_ADMIXTURE_M1_k2.pdf"),
#        width = 4, height = 3)

# M2
(admix_m2_ampse <- admix_ampse %>%
    filter(focus_group == "m2gh") %>%
    mutate(names = factor(names, level = order_m2)) %>%
    admix_fct(hide_legend = T))

# ggsave(plot = admix_m2_ampse, filename = paste0(path_figures, "MSAS_ADMIXTURE_M2_k2.pdf"),
#        width = 4, height = 3)

# M3
(admix_m3_ampse <- admix_ampse %>%
    filter(focus_group == "m3gh") %>%
    mutate(names = factor(names, level = order_m3)) %>%
    admix_fct(hide_legend = T))

# ggsave(plot = admix_m3_ampse, filename = paste0(path_figures, "MSAS_ADMIXTURE_M3_k2.pdf"),
#        width = 4, height = 3)


# ADMIXTURE MSAS Bootstrapping --------------------------------------------

# List Q_se files
listed_files_se_ampse <- list.files(path = paste0(path_admix_ampse, "bootstrapping/"),
                                    pattern = ".Q_se$",
                                    full.names = T)

# Read listed files and set names
files_se_ampse <- map(listed_files_se_ampse, read.table)
names(files_se_ampse) <- listed_files_se_ampse %>% 
  gsub(".*topweide_", "", .) %>% 
  gsub("\\.Q_se", "", .) %>% 
  gsub("\\.", "_", .)

# Adapt dfs in list
files_m1_se_ampse <- map(files_se_ampse[c(na.omit(str_extract(names(files_se_ampse),"m1_.*")))], 
                         function(df) df %>% 
                           mutate(names = str_extract(sample_names_admix_ampse$names, ".*(Ari|Rep).*") %>% 
                                    na.omit() %>% 
                                    c()))
files_m2_se_ampse <- map(files_se_ampse[c(na.omit(str_extract(names(files_se_ampse),"m2.*")))], 
                         function(df) df %>% 
                           mutate(names = str_extract(sample_names_admix_ampse$names, ".*(Ari|Art).*") %>% 
                                    na.omit() %>% 
                                    c()))
files_m3_se_ampse <- map(files_se_ampse[c(na.omit(str_extract(names(files_se_ampse),"m3.*")))], 
                         function(df) df %>% 
                           mutate(names = str_extract(sample_names_admix_ampse$names, ".*(Ari|Alg).*") %>% 
                                    na.omit() %>% 
                                    c()))

files_kcol_se_ampse <- c(files_m1_se_ampse, files_m2_se_ampse, files_m3_se_ampse)
files_kcol_se_ampse <- map2(files_kcol_se_ampse, names(files_kcol_se_ampse),
                            function(df, df_name) df %>%
                              mutate(n_clusters = as.numeric(str_extract(df_name, "(?<=_).*")),
                                     # add column with information about the samples included in the ADMIXTURE analysis
                                     focus_group = str_extract(df_name, ".*(?=_)")) %>%
                              left_join(sample_names_admix_ampse, by = "names") %>% 
                              pivot_longer(cols = contains("V"), names_to = "Cluster_name", values_to = "Cluster_value_SE") %>% 
                              mutate(Cluster_name = gsub("V", "Cluster_", Cluster_name) %>%
                                       factor(levels = paste0("Cluster_", 1:max(n_clusters)))))

# Bind all dfs and compute confidence interval (CI)
admix_se_ampse <- bind_rows(files_kcol_se_ampse) %>% 
  mutate(Cluster_value_CI = Cluster_value_SE*qnorm(0.975)) %>% 
  mutate(names = case_when(str_detect(names, "LT-Rep50-3") ~ "LT-Rep50-1",
                           str_detect(names, "Rep-3") ~ "Rep-2",
                           str_detect(names, "Art25-3") ~ "Art25-1",
                           TRUE ~ names))


# Join admix and admix_se
admix_all_ampse <- left_join(admix_ampse, admix_se_ampse, by = intersect(names(admix_ampse), names(admix_se_ampse)))

# M1
(admix_m1_ci_ampse <- admix_all_ampse %>% 
    filter(focus_group == "m1") %>% 
    mutate(names = factor(names, level = order_m1),
           Cluster_value_CI = case_when(Cluster_name == "Cluster_2" ~ Cluster_value_CI,
                                        TRUE ~ NaN)) %>%
    admix_fct(hide_legend = T) +
    geom_errorbar(aes(ymin = Cluster_value-Cluster_value_CI, ymax = Cluster_value+Cluster_value_CI), 
                  width = .15, linewidth = .3))

# ggsave(plot = admix_m1_ci_ampse, filename = paste0(path_figures, "Ampse_ADMIXTURE_M1_k2_CI.pdf"),
#        width = 4, height = 3)

# M2
(admix_m2_ci_ampse <- admix_all_ampse %>% 
    filter(focus_group == "m2gh") %>% 
    mutate(names = factor(names, level = order_m2),
           Cluster_value_CI = case_when(Cluster_name == "Cluster_2" ~ Cluster_value_CI,
                                        TRUE ~ NaN)) %>%
    admix_fct(hide_legend = T) +
    geom_errorbar(aes(ymin = Cluster_value-Cluster_value_CI, ymax = Cluster_value+Cluster_value_CI), 
                  width = .15, linewidth = .3))

# ggsave(plot = admix_m2_ci_ampse, filename = paste0(path_figures, "MSAS_ADMIXTURE_M2_k2_CI.pdf"),
#        width = 4, height = 3)

# M3
(admix_m3_ci_ampse <- admix_all_ampse %>% 
    filter(focus_group == "m3gh") %>% 
    mutate(names = factor(names, level = order_m3),
           Cluster_value_CI = case_when(Cluster_name == "Cluster_2" ~ Cluster_value_CI,
                                        TRUE ~ NaN)) %>%
    admix_fct(hide_legend = T) +
    geom_errorbar(aes(ymin = Cluster_value-Cluster_value_CI, ymax = Cluster_value+Cluster_value_CI), 
                  width = .15, linewidth = .3))

# ggsave(plot = admix_m3_ci_ampse, filename = paste0(path_figures, "MSAS_ADMIXTURE_M3_k2_CI.pdf"),
#        width = 4, height = 3)


# STRUCTURE MSAS ----------------------------------------------------------

# List files without priors
listed_files_strct_ampse <- list.files(path = path_strct_ampse,
                                       pattern = "_f$",
                                       full.names = T)

# Load and modify files without priors
files_strct_ampse <- map2(listed_files_strct_ampse, 
                          str_extract(listed_files_strct_ampse, "(?<=topweide_).*(?=_k)"),
                          function(file_name, df_name){
                            read.table(file_name, sep = "*") %>% # sep doesn't appear once -> reads files in one column
                              filter(str_detect(V1, "(P|M|F)-")) %>%
                              separate(V1, c("sample", "cluster"), sep = ":") %>% 
                              mutate(sample = str_extract(sample, "(P|M|F).*-\\d"),
                                     n_clusters = str_count(cluster, "\\d\\s\\d")+1,
                                     cluster = str_extract(cluster, "\\d.*\\d")) %>% 
                              separate(cluster, sep = " ", into = paste0("Cluster_", 1:.[1,"n_clusters"])) %>%
                              pivot_longer(cols = starts_with("Cluster"), names_to = "Cluster_name", values_to = "Cluster_value") %>% 
                              mutate(Cluster_value = as.numeric(Cluster_value),
                                     # Add focus group based on name of df
                                     focus_group = df_name)
                          })

# Load file with priors (manually created)
strct_priors_ampse <- read.table(paste0(path_strct_ampse, "ampse_topweide_priors_man.txt"),
                                 fill = T, header = T) %>% 
  pivot_longer(cols = starts_with("Cluster"), names_to = "Cluster_name", values_to = "Cluster_value") %>% 
  na.omit() %>% 
  # Add column with nr. of clusters
  group_by(focus_group) %>% 
  mutate(n_clusters = max(as.numeric(str_extract(Cluster_name, "\\d"))))



# Bind all dfs
strct_ampse <- bind_rows(files_strct_ampse, strct_priors_ampse) %>%
  left_join(samples_names_ampse, by = "sample")

# Plot
# M1
(strct_m1_ampse <- strct_ampse %>% 
    filter(focus_group == "m1ghf") %>% 
    mutate(names = factor(names, level = order_m1_str)) %>%
    admix_fct(hide_legend = T))

# ggsave(plot = strct_m1_ampse, filename = paste0(path_figures, "MSAS_STRUCTURE_M1_k2.pdf"),
#        width = 4, height = 3)

# M2
(strct_m2_ampse <- strct_ampse %>% 
    filter(focus_group == "m2gh") %>% 
    mutate(names = factor(names, level = order_m2)) %>%
    admix_fct(hide_legend = T))

# ggsave(plot = strct_m2_ampse, filename = paste0(path_figures, "MSAS_STRUCTURE_M2_k2.pdf"),
#        width = 4, height = 3)

# M3
(strct_m3_ampse <- strct_ampse %>% 
    filter(focus_group == "m3gh") %>% 
    mutate(names = factor(names, level = order_m3)) %>%
    admix_fct(hide_legend = T))

# ggsave(plot = strct_m3_ampse, filename = paste0(path_figures, "MSAS_STRUCTURE_M3_k2.pdf"),
#        width = 4, height = 3)
