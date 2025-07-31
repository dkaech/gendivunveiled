# Purpose: Load MSAS VCF file, create FST-heatmap and pairplot of allele frequency
# Author: Damian Käch <dkaech@ethz.ch>

# Load data ---------------------------------------------------------------


# Sample names based on vcf file
samples_names_ampse <- read.table(paste0(path_vcf, "MSAS/extended_loper/ampse_vcf_samples.txt"),
                                  col.names = "sample") %>% 
  mutate(POOLSIZE = case_when(str_detect(sample, "Ara|Ari|Rep") ~ 60,
                              str_detect(sample, "Art|Arc|Alg|M-1|F") ~ 120,
                              TRUE ~ 180),
         ploidy = case_when(str_detect(sample, "Ara|Ari|Rep|M-1|F") ~ 2,
                            str_detect(sample, "Art|Arc|Alg|") ~ 4,
                            TRUE ~ 6),
         names = case_when(str_detect(sample, "P") ~ str_remove(sample, "P-"),
                           str_detect(sample, "M-1") ~ str_replace(sample, "M-1-", "Rep"),
                           str_detect(sample, "F-LT") ~ str_replace(sample, "F-LT", "LT-Rep50"),
                           str_detect(sample, "F-EB") ~ str_replace(sample, "F-EB", "EB-Rep50"),
                           str_detect(sample, "M-2") ~ str_replace(sample, "M-2-", "Art"),
                           str_detect(sample, "M-3") ~ str_replace(sample, "M-3-", "Alg")) %>% 
           str_replace("75", "25"))

# Read VCF
ampse_raw <- vcf2pooldata(vcf.file = paste0(path_vcf, "MSAS/extended_loper/MSAS_loper.vcf"),
                          poolnames = samples_names_ampse$names,
                          poolsizes = samples_names_ampse$POOLSIZE)

# Subset pooldata file
ampse_pure <- pooldata.subset(ampse_raw, pool.index = which(str_detect(ampse_raw@poolnames, "^.{5}$")),
                              verbose = F)


# FST ---------------------------------------------------------------------



# Pairwise FST
ampse.pairwisefst_pure <- compute.pairwiseFST(ampse_pure)

# Max and min Fst
cat("Max FST:", max(ampse.pairwisefst_pure@PairwiseFSTmatrix, na.rm = T),
    "\nMin FST:", min(ampse.pairwisefst_pure@PairwiseFSTmatrix, na.rm = T),
    sep = " ")

# FST heatmap
(
  ampse_heatmap_fst <- ampse.pairwisefst_pure@PairwiseFSTmatrix %>% 
    # negative FST values are seen as zero values
    abs() %>% 
    .[order(rownames(.)),order(colnames(.))] %>% 
    get_lower_tri() %>% 
    reshape2::melt(na.rm = T) %>%
    mutate(Var1 = str_replace(Var1, "Rep-3", "Rep-2"),
           Var2 = str_replace(Var2, "Rep-3", "Rep-2")) %>% 
    ggplot(aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = value)) +
    scale_y_discrete(position = "right") +
    coord_flip() +
    scale_fill_gradient2(high = "#E69F00", mid = "#56B4E9",
                         limits = c(0,0.21), breaks = c(0,.05,.1,.15,.2),
                         name = bquote(italic("F")[ST] ~ "value")) +
    theme_minimal() +
    theme(text = element_text(color = "black"), 
          axis.text = element_text(size = 18),
          axis.text.x = element_text(angle = 70, hjust = .1),
          axis.title = element_blank(),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 13))
)

# ggsave(plot = ampse_heatmap_fst, filename = paste0(path_figures, "MSAS_FST_heatmap_pure.pdf"),
#        width = 11, height = 10)

# Save without legend
# ggsave(plot = ampse_heatmap_fst +
#          theme(legend.position = "none"),
#        filename = paste0(path_figures, "MSAS_FST_heatmap_pure_nole.pdf"),
#        width = 7, height = 7)


# Allele frequencies MSAS -------------------------------------------------

# Compute allele frequencies
ampse_frq <- ampse_raw@refallele.readcount/ampse_raw@readcoverage
colnames(ampse_frq) <- ampse_raw@poolnames
rownames(ampse_frq) <- paste(ampse_raw@snp.info$Chromosome, 
                             ampse_raw@snp.info$Position,
                             sep = "_")
ampse_frq <- na.omit(ampse_frq)
ampse_frq <- ampse_frq[,sort(colnames(ampse_frq))]
ampse_frq %>% head(3)
# Remove F-LT-1
ampse_frq <- ampse_frq[,colnames(ampse_frq) != "LT-Rep50-1"]

# Correlation allele frequencies
ampse_frq_df <- as.data.frame(ampse_frq) %>% 
  mutate(SNP = rownames(.)) %>% 
  pivot_longer(cols = -SNP, names_to = "Sample", values_to = "allele_freq") %>% 
  mutate(sample_grp = str_remove(Sample, "(-1|-2|-3)$"),
         sample_rep = str_extract(Sample, "(1|2|3)$")) %>% 
  select(-Sample) %>% 
  pivot_wider(names_from = sample_rep, names_prefix = "Replicate", values_from = allele_freq)

# Pairplot of biological replicates
ampse_frq_df %>% 
  filter(str_detect(.$sample_grp, "^.{3}$|^.{8}$")) %>% 
  .[,3:5] %>% 
  pairplot(save = F, custom_alpha = 0.1)

# ggsave(filename = paste0(path_figures, "MSAS_AF_biol.pdf"),
#        width = 16, height = 9, device = "pdf")

# Pairplot of technical replicates
ampse_frq_df %>% 
  filter(!str_detect(.$sample_grp, "^.{3}$|^.{8}$")) %>% 
  .[,3:5] %>% 
  pairplot(save = F, custom_alpha = 0.1)

# ggsave(filename = paste0(path_figures, "MSAS_AF_tech.pdf"),
#        width = 16, height = 9, device = "pdf")

