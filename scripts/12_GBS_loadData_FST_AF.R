# Purpose: Load GBS VCF file, create FST-heatmap and pairplot of allele frequency
# Author: Damian Käch <dkaech@ethz.ch>


# Load data ---------------------------------------------------------------

# Sample names based on vcf file
samples_names <- read.table(paste0(path_vcf, "GBS/sample_names_GBS.txt"),
                            col.names = "sample") %>% 
  mutate(poolsize = case_when(str_detect(sample, "Ara|Ari|Rep") ~ 60,
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
gbs_raw <- vcf2pooldata(vcf.file = paste0(path_vcf, "GBS/gbs.vcf"),
                        poolnames = samples_names$names,
                        poolsizes = samples_names$poolsize)

# Subset pooldata file
gbs_pure <- pooldata.subset(gbs_raw, pool.index = which(str_detect(gbs_raw@poolnames, "^.{5}$")),
                            verbose = F)


# FST ---------------------------------------------------------------------


# Pairwise FST
gbs.pairwisefst_pure <- compute.pairwiseFST(gbs_pure)



# Max and min Fst
cat("Max FST:", max(gbs.pairwisefst_pure@PairwiseFSTmatrix, na.rm = T),
    "\nMin FST:", min(gbs.pairwisefst_pure@PairwiseFSTmatrix, na.rm = T),
    sep = " ")

# FST heatmap
(
  gbs_heatmap_fst <- gbs.pairwisefst_pure@PairwiseFSTmatrix %>% 
    .[order(rownames(.)),order(colnames(.))] %>% 
    get_lower_tri() %>% 
    reshape2::melt(na.rm = T) %>%
    ggplot(aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = value)) +
    scale_y_discrete(position = "right") +
    coord_flip() +
    scale_fill_gradient2(high = "#E69F00", mid = "#56B4E9",
                         limits = c(0,.21), breaks = c(0,.05,.1,.15,.2),
                         name = bquote(italic("F")[ST] ~ "value")) +
    theme_minimal() +
    theme(text = element_text(color = "black"), 
          axis.text = element_text(size = 18),
          axis.text.x = element_text(angle = 70, hjust = .1),
          axis.title = element_blank(),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 13))
)

# ggsave(plot = gbs_heatmap_fst, filename = paste0(path_figures, "GBS_FST_heatmap_pure.pdf"),
#        width = 11, height = 10)

# Save without legend
# ggsave(plot = gbs_heatmap_fst +
#          theme(legend.position = "none"),
#        filename = paste0(path_figures, "GBS_FST_heatmap_pure_nole.pdf"),
#        width = 7, height = 7)

# Allele frequencies ------------------------------------------------------

# Compute allele frequencies
gbs_frq <- gbs_raw@refallele.readcount/gbs_raw@readcoverage
colnames(gbs_frq) <- gbs_raw@poolnames
rownames(gbs_frq) <- paste(gbs_raw@snp.info$Chromosome, 
                           gbs_raw@snp.info$Position,
                           sep = "_")

gbs_frq <- na.omit(gbs_frq)
gbs_frq <- gbs_frq[,sort(colnames(gbs_frq))]
gbs_frq %>% head(3)
# Remove F-LT-1
gbs_frq <- gbs_frq[,colnames(gbs_frq) != "LT-Rep50-1"]

# Correlation allele frequencies
gbs_frq_df <- as.data.frame(gbs_frq) %>% 
  mutate(SNP = rownames(.)) %>% 
  pivot_longer(cols = -SNP, names_to = "Sample", values_to = "allele_freq") %>% 
  mutate(sample_grp = str_remove(Sample, "(-1|-2|-3)$"),
         sample_rep = str_extract(Sample, "(1|2|3)$")) %>% 
  select(-Sample) %>% 
  pivot_wider(names_from = sample_rep, names_prefix = "Replicate", values_from = allele_freq)

# Pairplot of biological replicates (takes a while)
gbs_frq_df %>% 
  filter(str_detect(.$sample_grp, "^.{3}$|^.{8}$")) %>% 
  # Select only numerical columns
  .[,3:5] %>% 
  pairplot(save = F, custom_alpha = 0.002)

# Save as pdf (45MB)
# ggsave(filename = paste0(path_figures, "GBS_AF_biol.pdf"),
#        width = 16, height = 9, device = "pdf")
# Save as jpg
# ggsave(filename = paste0(path_figures, "GBS_AF_biol.jpg"),
#        width = 16, height = 9, device = "jpg", dpi = 100)

# Pairplot of technical replicates (takes a while)
gbs_frq_df %>% 
  filter(!str_detect(.$sample_grp, "^.{3}$|^.{8}$")) %>% 
  # Select only numerical columns
  .[,3:5] %>% 
  pairplot(save = F, custom_alpha = 0.002)

# Save as pdf (38 MB)
# ggsave(filename = paste0(path_figures, "GBS_AF_tech.pdf"),
#        width = 16, height = 9, device = "pdf")
# Save as jpg
# ggsave(filename = paste0(path_figures, "GBS_AF_tech.jpg"),
#        width = 16, height = 9, device = "jpg", dpi = 100)