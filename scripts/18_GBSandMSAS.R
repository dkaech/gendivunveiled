# Purpose: Analyse GBS and MSAS combined
# Author: Damian Käch <dkaech@ethz.ch>


# Extract legends ---------------------------------------------------------

# Legends are extracted to add them along the no-legend figures

# Extract legend of FST heatmap
legend_fst_heatmap <- get_legend(gbs_heatmap_fst)

# ggsave(plot = legend_fst_heatmap, filename = paste0(path_figures, "FST_heatmap_legend.pdf"),
#        width = 1, height = 1.7)

# Extract legend of ADMIXTURE
admix_legend <- admix %>% 
  filter(focus_group == "m1") %>% 
  mutate(names = factor(names, level = order_m1),
         Cluster_name = gsub("_", " ", Cluster_name)) %>%
  admix_fct() +
  theme(legend.position = "bottom",
        legend.title = element_blank())
admix_legend <- get_legend(admix_legend)

# ggsave(plot = admix_legend, filename = paste0(path_figures, "GBS_ADMIXTURE_legend.pdf"),
#        width = 2, height = .5)


# ADMIXTURE CV ------------------------------------------------------------

# Comparison of cross-validation error between GBS and MSAS (for mixtures)

listed_files_cv <- list.files(path = paste0(path_admix, "bootstrapping"),
                              pattern = "*error",
                              full.names = T)

files_cv <- map(listed_files_cv, read.table)
names(files_cv) <- listed_files_cv %>% 
  gsub(".*topweide_", "", .) %>% 
  gsub("\\.cv\\.error", "", .)

# Adapt dfs in list
df_cv <- map2(files_cv, names(files_cv),
              function(df, df_name) df %>%
                rename(n_clusters = V1,
                       cluster_error = V2) %>% 
                mutate(Mixture = str_extract(df_name, "[:digit:]") %>% 
                         paste0("M", .),
                       approach = "GBS")) %>% 
  bind_rows()

listed_files_cv_ampse <- list.files(path = paste0(path_admix_ampse, "bootstrapping"),
                                    pattern = "*error",
                                    full.names = T)

files_cv_ampse <- map(listed_files_cv_ampse, read.table)
names(files_cv_ampse) <- listed_files_cv_ampse %>% 
  gsub(".*topweide_", "", .) %>% 
  gsub("\\.cv\\.error", "", .)

# Adapt dfs in list
df_cv_ampse <- map2(files_cv_ampse, names(files_cv_ampse),
                    function(df, df_name) df %>%
                      rename(n_clusters = V1,
                             cluster_error = V2) %>% 
                      mutate(Mixture = str_extract(df_name, "[:digit:]") %>% 
                               paste0("M", .),
                             approach = "MSAS")) %>% 
  bind_rows()

df_cv_all <- bind_rows(df_cv, df_cv_ampse)

(admix_cv <- ggplot(df_cv_all, aes(x = n_clusters, y = cluster_error, color = Mixture)) +
    geom_line(linewidth = .5) +
    geom_point() +
    scale_color_manual(values = c("#E69F00","#56B4E9","#D55E00")) +
    facet_grid(~approach) +
    labs(x = "K", y = "Cross-validation error") +
    theme_bw() +
    theme(text = element_text(colour = "black", size = 11),
          strip.text = element_text(color = "black", size = 12)))

# ggsave(plot = admix_cv, filename = paste0(path_figures, "ADMIXTURE_CV.pdf"),
#        width = 6, height = 3)


# ADMIXTURE K>2 -----------------------------------------------------------

# Comparison of cluster memberships between GBS and MSAS (for pure samples)

palette_admix <- hcl.colors(max(admix$n_clusters), "Roma")

# Plot one of pure, m1, m2, m3
admix_focus_fct <- function(df){
  df %>% 
    mutate(n_clusters = paste0("K = ", n_clusters)) %>% 
    ggplot(mapping = aes(x = names, y = Cluster_value, fill = Cluster_name)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = palette_admix) +
    facet_grid(n_clusters ~ approach) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(color = "black", size = 11),
          strip.text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.title = element_blank())
}

admix_gbs_ampse <- admix %>% 
  mutate(approach = "GBS") %>% 
  bind_rows(admix_ampse %>% 
              mutate(approach = "MSAS")) %>% 
  mutate(focus_group = case_when(str_detect(focus_group, "p") ~ "pure",
                                 TRUE ~ focus_group))

# Print and save plots
# Pure
admix_focus_fct(admix_gbs_ampse %>% 
                  mutate(names = factor(names, level = c(paste0("Ara-", 1:3),
                                                         paste0("Ari-", 1:3),
                                                         paste0("Rep-", 1:3),
                                                         paste0("Art-", 1:3),
                                                         paste0("Arc-", 1:3),
                                                         paste0("Alg-", 1:3)))) %>% 
                  filter(focus_group == "pure"))
# ggsave(filename = paste0(path_figures, "ADMIXTURE_pure_k_all.pdf"),
#        width = 8, height = 10, device = "pdf")
# M1
admix_focus_fct(admix_gbs_ampse %>% 
                  mutate(names = factor(names, level = order_m1)) %>% 
                  filter(focus_group == "m1"))
# ggsave(filename = paste0(path_figures, "ADMIXTURE_m1_k_all.pdf"),
#        width = 8, height = 10, device = "pdf")
# M2
admix_focus_fct(admix_gbs_ampse %>% 
                  mutate(names = factor(names, level = order_m2)) %>% 
                  filter(str_detect(focus_group, "m2")))
# ggsave(filename = paste0(path_figures, "ADMIXTURE_m2_k_all.pdf"),
#        width = 8, height = 6, device = "pdf")
# M3
admix_focus_fct(admix_gbs_ampse %>% 
                  mutate(names = factor(names, level = order_m3)) %>% 
                  filter(str_detect(focus_group, "m3")))
# ggsave(filename = paste0(path_figures, "ADMIXTURE_m3_k_all.pdf"),
#        width = 8, height = 6, device = "pdf")

# STRUCTURE ---------------------------------------------------------------

# Comparison of cluster memberships between GBS and MSAS (for pure samples)

strct_gbs_ampse <- strct %>% 
  mutate(approach = "GBS") %>% 
  bind_rows(strct_ampse %>% 
              mutate(approach = "MSAS")) %>% 
  mutate(focus_group = case_when(str_detect(focus_group, "p") ~ "pure",
                                 TRUE ~ focus_group))

# Print and save plot with pure samples
admix_focus_fct(strct_gbs_ampse %>% 
                  mutate(names = factor(names, level = c(paste0("Ara-", 1:3),
                                                         paste0("Ari-", 1:3),
                                                         paste0("Rep-", 1:3),
                                                         paste0("Art-", 1:3),
                                                         paste0("Arc-", 1:3),
                                                         paste0("Alg-", 1:3)))) %>% 
                  filter(focus_group == "pure"))
# ggsave(filename = paste0(path_figures, "STRUCTURE_pure_k_all.pdf"),
#        width = 8, height = 10, device = "pdf")