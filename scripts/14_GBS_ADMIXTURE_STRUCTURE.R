# Purpose: Analyse ADMIXTURE and STRUCTURE output of GBS
# Author: Damian Käch <dkaech@ethz.ch>

# ADMIXTURE ---------------------------------------------------------------

# Read file with sample names (to get correct order of sample names)
sample_names_admix <- read.table(paste0(path_admix, "gbs_topweide_nolt.list"),
                                 sep = "\t") %>%
  rename(Sample = V1,
         Population = V2) %>%
  mutate(Sample = str_extract(Sample, "(P-|M-|F-).*(?=_R1)")) %>%
  left_join(select(samples_names, sample, names), by = c("Sample" = "sample"))

# List Q files
listed_files_focus <- list.files(path = paste0(path_admix, "focus_groups/"),
                                 pattern = ".Q$",
                                 full.names = T)

# Read listed files and set names
files_focus <- map(listed_files_focus, read.table)
names(files_focus) <- listed_files_focus %>%
  gsub(".*topweide_", "", .) %>%
  gsub("\\.Q", "", .) %>%
  gsub("\\.", "_", .)

# Adapt dfs in list
files_pure <- map(files_focus[c(na.omit(str_extract(names(files_focus),"ghpure.*")))], 
                  function(df) df %>% 
                    mutate(names = str_extract(sample_names_admix$names, "^...-.*") %>% 
                             na.omit() %>% 
                             c()))
files_m1 <- map(files_focus[c(na.omit(str_extract(names(files_focus),"m1_.*")))],
                function(df) df %>%
                  mutate(names = str_extract(sample_names_admix$names, ".*(Ari|Rep).*") %>%
                           na.omit() %>%
                           c()))
files_m2 <- map(files_focus[c(na.omit(str_extract(names(files_focus),"m2.*")))],
                function(df) df %>%
                  mutate(names = str_extract(sample_names_admix$names, ".*(Ari|Art).*") %>%
                           na.omit() %>%
                           c()))
files_m3 <- map(files_focus[c(na.omit(str_extract(names(files_focus),"m3.*")))],
                function(df) df %>%
                  mutate(names = str_extract(sample_names_admix$names, ".*(Ari|Alg).*") %>%
                           na.omit() %>%
                           c()))

files_kcol <- c(files_pure, files_m1, files_m2, files_m3)
files_kcol <- map2(files_kcol, names(files_kcol),
                   function(df, df_name) df %>%
                     mutate(n_clusters = as.numeric(str_extract(df_name, "(?<=_).*")),
                            # add column with information about the samples included in the ADMIXTURE analysis
                            focus_group = str_extract(df_name, ".*(?=_)")) %>%
                     left_join(sample_names_admix, by = "names") %>%
                     pivot_longer(cols = contains("V"), names_to = "Cluster_name", values_to = "Cluster_value") %>%
                     mutate(Cluster_name = gsub("V", "Cluster_", Cluster_name) %>%
                              factor(levels = paste0("Cluster_", 1:max(n_clusters)))))

# Bind all dfs
admix <- bind_rows(files_kcol) %>%
  mutate(names = str_replace(names, "LT-Rep50-3", "LT-Rep50-1"))

admix_fct <- function(df, hide_legend = FALSE){
  df %>%
    filter(n_clusters == 2) %>%
    ggplot(mapping = aes(x = names, y = Cluster_value, fill = Cluster_name)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = c("#E69F00","#56B4E9")) +
    scale_y_continuous(labels = scales::percent, breaks = c(0,.25,.5,.75,1)) +
    coord_cartesian(ylim = c(-0.01,1.01), expand = F) +
    theme_bw() +
    {if (hide_legend)theme(legend.position = "none")} +
    theme(text = element_text(color = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, size = 13),
          axis.text.y = element_text(size = 11),
          axis.title = element_blank(),
          axis.ticks.x = element_blank())
}

# M1
order_m1 <- c(paste0("Rep50-", 1:3), paste0("Rep25-", 1:3),
              paste0("EB-Rep50-", 1:3), paste0("LT-Rep50-", 1:2),
              paste0("Rep-", 1:3), paste0("Ari-", 1:3))

(admix_m1 <- admix %>%
    filter(focus_group == "m1") %>%
    mutate(names = factor(names, level = order_m1)) %>%
    admix_fct(hide_legend = T))

# ggsave(plot = admix_m1, filename = paste0(path_figures, "GBS_ADMIXTURE_M1_k2.pdf"),
#        width = 4, height = 3)

# M2
order_m2 <- c(paste0("Art50-", 1:3), paste0("Art25-", 1:3),
              paste0("Art-", 1:3), paste0("Ari-", 1:3))
(admix_m2 <- admix %>%
    filter(focus_group == "m2") %>%
    mutate(names = factor(names, level = order_m2),
           # Invert order of cluster names to have common cluster colors like M1 and M2
           Cluster_name = factor(Cluster_name, level = rev(levels(Cluster_name)))) %>%
    admix_fct(hide_legend = T))

# ggsave(plot = admix_m2, filename = paste0(path_figures, "GBS_ADMIXTURE_M2_k2.pdf"),
#        width = 4, height = 3)

# M3
order_m3 <- c(paste0("Alg50-", 1:3), paste0("Alg25-", 1:3),
              paste0("Alg-", 1:3), paste0("Ari-", 1:3))
(admix_m3 <- admix %>%
    filter(focus_group == "m3") %>%
    mutate(names = factor(names, level = order_m3)) %>%
    admix_fct(hide_legend = T))

# ggsave(plot = admix_m3, filename = paste0(path_figures, "GBS_ADMIXTURE_M3_k2.pdf"),
#        width = 4, height = 3)

# ADMIXTURE Bootstrapping -------------------------------------------------

# List Q_se files
listed_files_se <- list.files(path = paste0(path_admix, "bootstrapping/"),
                              pattern = ".Q_se$",
                              full.names = T)

# Read listed files and set names
files_se <- map(listed_files_se, read.table)
names(files_se) <- listed_files_se %>% 
  gsub(".*topweide_", "", .) %>% 
  gsub("\\.Q_se", "", .) %>% 
  gsub("\\.", "_", .)

# Adapt dfs in list
files_m1_se <- map(files_se[c(na.omit(str_extract(names(files_se),"m1_.*")))], 
                   function(df) df %>% 
                     mutate(names = str_extract(sample_names_admix$names, ".*(Ari|Rep).*") %>% 
                              na.omit() %>% 
                              c()))
files_m2_se <- map(files_se[c(na.omit(str_extract(names(files_se),"m2.*")))], 
                   function(df) df %>% 
                     mutate(names = str_extract(sample_names_admix$names, ".*(Ari|Art).*") %>% 
                              na.omit() %>% 
                              c()))
files_m3_se <- map(files_se[c(na.omit(str_extract(names(files_se),"m3.*")))], 
                   function(df) df %>% 
                     mutate(names = str_extract(sample_names_admix$names, ".*(Ari|Alg).*") %>% 
                              na.omit() %>% 
                              c()))

files_kcol_se <- c(files_m1_se, files_m2_se, files_m3_se)
files_kcol_se <- map2(files_kcol_se, names(files_kcol_se),
                      function(df, df_name) df %>%
                        mutate(n_clusters = as.numeric(str_extract(df_name, "(?<=_).*")),
                               # add column with information about the samples included in the ADMIXTURE analysis
                               focus_group = str_extract(df_name, ".*(?=_)")) %>%
                        left_join(sample_names_admix, by = "names") %>% 
                        pivot_longer(cols = contains("V"), names_to = "Cluster_name", values_to = "Cluster_value_SE") %>% 
                        mutate(Cluster_name = gsub("V", "Cluster_", Cluster_name) %>%
                                 factor(levels = paste0("Cluster_", 1:max(n_clusters)))))

# Bind all dfs and compute confidence interval (CI)
admix_se <- bind_rows(files_kcol_se) %>% 
  mutate(Cluster_value_CI = Cluster_value_SE*qnorm(0.975)) %>% 
  mutate(names = str_replace(names, "LT-Rep50-3", "LT-Rep50-1"))


# Join admix and admix_se
admix_all <- left_join(admix, admix_se, by = intersect(names(admix), names(admix_se)))

# M1
(admix_m1_ci <- admix_all %>% 
    filter(focus_group == "m1") %>% 
    mutate(names = factor(names, level = order_m1),
           Cluster_value_CI = case_when(Cluster_name == "Cluster_2" ~ Cluster_value_CI,
                                        TRUE ~ NaN)) %>%
    admix_fct(hide_legend = T) +
    geom_errorbar(aes(ymin = Cluster_value-Cluster_value_CI, ymax = Cluster_value+Cluster_value_CI), 
                  width = .15, linewidth = .3))

# ggsave(plot = admix_m1_ci, filename = paste0(path_figures, "GBS_ADMIXTURE_M1_k2_CI.pdf"),
#        width = 4, height = 3)

# M2
(admix_m2_ci <- admix_all %>% 
    filter(focus_group == "m2") %>% 
    mutate(names = factor(names, level = order_m2),
           Cluster_value_CI = case_when(Cluster_name == "Cluster_2" ~ Cluster_value_CI,
                                        TRUE ~ NaN),
           # Ari is in Cluster 2 -> switch clusters to be consistent with other results
           Cluster_name = case_when(str_detect(Cluster_name, "Cluster_1") ~ "Cluster_99",
                                    str_detect(Cluster_name, "Cluster_2") ~ "Cluster_1",
                                    str_detect(Cluster_name, "Cluster_99") ~ "Cluster_2")) %>%
    admix_fct(hide_legend = T) +
    # Ari is in Cluster 2 -> compute 1-... to correctly place CIs
    geom_errorbar(aes(ymin = 1-Cluster_value-Cluster_value_CI, ymax = 1-Cluster_value+Cluster_value_CI), 
                  width = .15, linewidth = .3))

# ggsave(plot = admix_m2_ci, filename = paste0(path_figures, "GBS_ADMIXTURE_M2_k2_CI.pdf"),
#        width = 4, height = 3)

# M3
(admix_m3_ci <- admix_all %>% 
    filter(focus_group == "m3") %>% 
    mutate(names = factor(names, level = order_m3),
           Cluster_value_CI = case_when(Cluster_name == "Cluster_2" ~ Cluster_value_CI,
                                        TRUE ~ NaN)) %>%
    admix_fct(hide_legend = T) +
    geom_errorbar(aes(ymin = Cluster_value-Cluster_value_CI, ymax = Cluster_value+Cluster_value_CI), 
                  width = .15, linewidth = .3))

# ggsave(plot = admix_m3_ci, filename = paste0(path_figures, "GBS_ADMIXTURE_M3_k2_CI.pdf"),
#        width = 4, height = 3)


# STRUCTURE ---------------------------------------------------------------

# List files without priors
listed_files_strct <- list.files(path = path_strct,
                                 pattern = "_f$",
                                 full.names = T)

# Load and modify files without priors
files_strct <- map2(listed_files_strct,
                    str_extract(listed_files_strct, "(?<=topweide_).*(?=_k)"),
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
strct_priors <- read.table(paste0(path_strct, "gbs_topweide_priors_man.txt"),
                           fill = T, header = T) %>% 
  pivot_longer(cols = starts_with("Cluster"), names_to = "Cluster_name", values_to = "Cluster_value") %>% 
  na.omit() %>% 
  # Add column with nr. of clusters
  group_by(focus_group) %>% 
  mutate(n_clusters = max(as.numeric(str_extract(Cluster_name, "\\d"))))



# Bind all dfs
strct <- bind_rows(files_strct, strct_priors) %>%
  left_join(sample_names_admix, by = c("sample" = "Sample"))


# Plot
# M1
order_m1_str <- c(paste0("Rep50-", 1:3), paste0("Rep25-", 1:3),
                  paste0("EB-Rep50-", 1:3), paste0("LT-Rep50-", 2:3),
                  paste0("Rep-", 1:3), paste0("Ari-", 1:3))

(strct_m1 <- strct %>% 
    filter(focus_group == "m1ghf") %>% 
    mutate(names = factor(names, level = order_m1_str)) %>%
    admix_fct(hide_legend = T))

# ggsave(plot = strct_m1, filename = paste0(path_figures, "GBS_STRUCTURE_M1_k2.pdf"),
#        width = 4, height = 3)

# M2
(strct_m2 <- strct %>% 
    filter(focus_group == "m2gh") %>% 
    mutate(names = factor(names, level = order_m2)) %>%
    admix_fct(hide_legend = T))

# ggsave(plot = strct_m2, filename = paste0(path_figures, "GBS_STRUCTURE_M2_k2.pdf"),
#        width = 4, height = 3)

# M3
(strct_m3 <- strct %>% 
    filter(focus_group == "m3gh") %>% 
    mutate(names = factor(names, level = order_m3)) %>%
    admix_fct(hide_legend = T))

# ggsave(plot = strct_m3, filename = paste0(path_figures, "GBS_STRUCTURE_M3_k2.pdf"),
#        width = 4, height = 3)

