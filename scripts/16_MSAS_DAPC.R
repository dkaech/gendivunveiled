# Purpose: Analyse allele frequency of MSAS using discriminant analysis of 
#           principal components (DAPC)
# Author: Damian Käch <dkaech@ethz.ch>


# DAPC --------------------------------------------------------------------


# Prepare data

ampse_frq_dapc <- t(ampse_frq) %>% 
  as.data.frame() %>% 
  mutate(Sample = rownames(.),
         Group = gsub("-.$", "", Sample) %>% 
           factor(),
         Group_nr = as.numeric(Group)) %>% 
  select(Sample, Group, Group_nr, contains("_"))

ampse_frq_dapc_pure <- ampse_frq_dapc %>% 
  filter(str_detect(Sample, "^.{5}$")) %>% 
  mutate(Group = as.factor(as.character(Group)))

ampse_frq_dapc_mix <- ampse_frq_dapc %>% 
  filter(!str_detect(Sample, "^.{5}$")) %>% 
  mutate(Group = as.factor(as.character(Group)) )%>% 
  select(-c("Sample","Group","Group_nr"))

# Reproducibility
set.seed(11)

# Training and cross-validaton (CV) based on pure samples

# The following line doesn't need to be run. User can use the load() function
# below as the object was saved before.
# system.time(
#   pramx_pure_ampse <- xvalDapc(ampse_frq_dapc_pure %>%
#                            # Keep only numerical variables
#                            select(-c("Sample","Group","Group_nr")),
#                          grp = ampse_frq_dapc_pure$Group,
#                          training.set = .9,
#                          result = "groupMean",
#                          # All variables have the same range of 0 to 1
#                          # --> no centering and no scaling
#                          center = F,
#                          scale = F,
#                          n.pca = NULL,
#                          n.pca.max = 18,
#                          n.rep = 30,
#                          xval.plot = T)
# )
# ~2 seconds

# Save and/or load the cross-validation object
# save(pramx_pure_ampse, file = paste0(path_tables, "MSAS_dapc_CV_pure.RData"))
load(file = paste0(path_tables, "MSAS_dapc_CV_pure.RData"))

pramx_pure_ampse$`Number of PCs Achieving Lowest MSE`
# "4"

# Predict mixtures based on training
pred.mix_ampse <- predict.dapc(pramx_pure_ampse$DAPC,
                               newdata = ampse_frq_dapc_mix)

df_grp2_1_ampse <- pramx_pure_ampse$DAPC$grp.coord %>% 
  as.data.frame() %>% 
  mutate(Name = rownames(.),
         Category = "Group")

df_grp2_2_ampse <- pred.mix_ampse$ind.scores %>% 
  as.data.frame() %>% 
  mutate(Name = rownames(.), 
         Name = gsub("-1$|-2$|-3$", "", Name)) %>% 
  group_by(Name) %>% 
  summarise(across(starts_with("LD"), mean)) %>% 
  mutate(Category = "Group") %>% 
  ungroup()

df_grp2_ampse <- bind_rows(df_grp2_1_ampse, df_grp2_2_ampse) %>% 
  mutate(Colors_3 = case_when(str_detect(Name, "^.{3}$") ~ "P",
                              str_detect(Name, "^.{5}$") ~ "M",
                              TRUE ~ "F"))

df_ind2_ampse <- pramx_pure_ampse$DAPC$ind.coord %>% 
  as.data.frame() %>% 
  bind_rows(as.data.frame(pred.mix_ampse$ind.scores)) %>% 
  mutate(Name = rownames(.),
         Category = case_when(str_detect(Name, "^.{5}$") ~ "Individual_pure",
                              TRUE ~ "Individual_mix"))

df_ind_grp2_ampse <- bind_rows(df_grp2_1_ampse, df_grp2_2_ampse, df_ind2_ampse) %>% 
  mutate(Group = gsub("-1$|-2$|-3$", "", Name),
         Colors_3 = case_when(str_detect(Name, "^.{3}$|^.{3}-.$") ~ "P",
                              str_detect(Name, "EB|LT") ~ "F",
                              TRUE ~ "M")
  )


dapc_ggplot_ampse <- ggplot(data = df_ind_grp2_ampse, aes(x = LD1, y = LD2)) +
  geom_point(aes(shape = Category, 
                 color = Colors_3, 
                 alpha = Category,
                 size = Category)) +
  scale_shape_manual(values = c("Group" = 4,
                                "Individual_pure" = 19,
                                "Individual_mix" = 15)) +
  scale_alpha_manual(values = c("Group" = 1,
                                "Individual_pure" = 0.3,
                                "Individual_mix" = 0.3)) +
  scale_size_manual(values = c("Group" = 5,
                               "Individual_pure" = 10,
                               "Individual_mix" = 10)) +
  geom_text_repel(data = df_grp2_ampse, 
                  aes(label = Name, color = Colors_3),
                  max.overlaps = 50, 
                  # Set distance between labels and points
                  box.padding = 1,
                  # Bellow following segment length, no
                  # segments are drawn (we don't want segments)
                  min.segment.length = 3,
                  size = 5.3) +
  scale_color_manual(values = c("F" = "#E69F00",
                                "P" = "#56B4E9",
                                "M" = "#D55E00")) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color = NA)
  )

dapc_ggplot_ampse

# ggsave(plot = dapc_ggplot_ampse,
#        filename = paste0(path_figures, "MSAS_DAPC_predmix_CV10.pdf"),
#        height = 7, width = 7, device = "pdf")

# Plot LD1 and LD2 for each mixture

# LD1
df_ari_ampse <- df_ind_grp2_ampse %>% 
  filter(str_detect(Name, "Ari")) %>% 
  mutate(focus_group = "M1")

df_ind_focus_ampse <- df_ind_grp2_ampse %>% 
  mutate(focus_group = case_when(str_detect(Name, "Rep") ~ "M1",
                                 str_detect(Name, "Art") ~ "M2",
                                 str_detect(Name, "Alg") ~ "M3",
                                 TRUE ~ "")) %>%
  filter(str_detect(focus_group, "M.*")) %>% 
  bind_rows(df_ari_ampse) %>% 
  bind_rows(df_ari_ampse %>% 
              mutate(focus_group = "M2")) %>% 
  bind_rows(df_ari_ampse %>% 
              mutate(focus_group = "M3"))

df_ari_txt_ampse <- df_grp2_ampse %>% 
  filter(str_detect(Name, "Ari$")) %>% 
  mutate(focus_group = "M1")

df_grp_txt_ampse <- df_grp2_ampse %>% 
  mutate(focus_group = case_when(str_detect(Name, "Rep") ~ "M1",
                                 str_detect(Name, "Art") ~ "M2",
                                 str_detect(Name, "Alg") ~ "M3",
                                 TRUE ~ "")) %>% 
  filter(str_detect(focus_group, "M")) %>% 
  bind_rows(df_ari_txt_ampse) %>% 
  bind_rows(df_ari_txt_ampse %>% 
              mutate(focus_group = "M2")) %>% 
  bind_rows(df_ari_txt_ampse %>% 
              mutate(focus_group = "M3"))

dapc_ggplot_ld_ampse <- ggplot(data = df_ind_focus_ampse,
                               aes(y = LD1, x = focus_group)) +
  geom_point(aes(shape = Category, 
                 color = Colors_3, 
                 alpha = Category,
                 size = Category)) +
  scale_shape_manual(values = c("Group" = 4,
                                "Individual_pure" = 19,
                                "Individual_mix" = 15)) +
  scale_alpha_manual(values = c("Group" = 1,
                                "Individual_pure" = 0.3,
                                "Individual_mix" = 0.3)) +
  scale_size_manual(values = c("Group" = 5,
                               "Individual_pure" = 10,
                               "Individual_mix" = 10)) +
  geom_text_repel(data = df_grp_txt_ampse, 
                  aes(label = Name, color = Colors_3),
                  max.overlaps = 50, 
                  # Set distance between labels and points
                  box.padding = 1,
                  # Bellow following segment length, no
                  # segments are drawn (we don't want segments)
                  min.segment.length = 2,
                  size = 5.3) +
  scale_color_manual(values = c("F" = "#E69F00",
                                "P" = "#56B4E9",
                                "M" = "#D55E00")) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(x = "Mixture") +
  # Same y limits like dapc_ggplot_ampse
  ylim(layer_scales(dapc_ggplot_ampse)$x$range$range) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color = NA)
  )


dapc_ggplot_ld_ampse

# ggsave(plot = dapc_ggplot_ld_ampse,
#        filename = paste0(path_figures, "MSAS_DAPC_predmix_CV10_LD1.pdf"),
#        height = 3, width = 7, device = "pdf")

# LD2
dapc_ggplot_ld2_ampse <- ggplot(data = df_ind_focus_ampse,
                                aes(y = LD2, x = focus_group)) +
  geom_point(aes(shape = Category, 
                 color = Colors_3, 
                 alpha = Category,
                 size = Category)) +
  scale_shape_manual(values = c("Group" = 4,
                                "Individual_pure" = 19,
                                "Individual_mix" = 15)) +
  scale_alpha_manual(values = c("Group" = 1,
                                "Individual_pure" = 0.3,
                                "Individual_mix" = 0.3)) +
  scale_size_manual(values = c("Group" = 5,
                               "Individual_pure" = 10,
                               "Individual_mix" = 10)) +
  geom_text_repel(data = df_grp_txt_ampse,
                  aes(label = Name, color = Colors_3),
                  max.overlaps = 50,
                  # Set distance between labels and points
                  box.padding = .8,
                  # Bellow following segment length, no
                  # segments are drawn (we don't want segments)
                  min.segment.length = 2,
                  size = 5.3) +
  scale_color_manual(values = c("F" = "#E69F00",
                                "P" = "#56B4E9",
                                "M" = "#D55E00")) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(x = "Mixture") +
  # Same y limits like dapc_ggplot_ampse
  ylim(layer_scales(dapc_ggplot_ampse)$y$range$range) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color = NA)
  )
dapc_ggplot_ld2_ampse

# ggsave(plot = dapc_ggplot_ld2_ampse,
#        filename = paste0(path_figures, "MSAS_DAPC_predmix_CV10_LD2.pdf"),
#        height = 7, width = 3, device = "pdf")
