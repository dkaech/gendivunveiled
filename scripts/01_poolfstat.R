# MSAS Poolfstat analysis with sp-sorted reads
# Author: Miguel Loera
# 29.05.2023

library(poolfstat)
library(tidyverse)
library(factoextra)
library(gridExtra)
library(ggforce)
library(reshape2)
library(broom)
library(AICcmodavg)
library(gplots)
library(dendextend)
library(ComplexHeatmap)
library(ggpubr)
library(ggridges)



# sample names ------------------------------------------------------------

results.final <- data.frame()
bxplots <- list()
fstplot <- list()

species.df <- data.frame(species = c('Dg', 'Fp', 'Lp', 'Tp', 'Tr'),
                         fspecies = c(
                           '*Dactylis glomerata* (MSAS)', 
                           '*Festuca pratensis* (MSAS)', 
                           '*Lolium perenne* (MSAS)', 
                           '*Trifol pratense* (MSAS)', 
                           '*Trifolium repens* (MSAS)'))



pool_size = data.frame(species =  c('Dg', 'Fp', 'Lp', 'Tp', 'Tr'), pool_size = c(80,40,40,40,80))


for (species in species.df$species) {

results <- data.frame(species = c(species))


samples.as <- data.frame(Sample_Name = c(paste0(sprintf("%s-A-", species), 1:3), paste0(sprintf("%s-B-", species), 1:3)))

samples.as <- samples.as %>% 
  mutate(species = str_extract(Sample_Name, "^[^-]+"),
         cultivar = str_extract(Sample_Name, "(?<=-)[^-]+"),
         rep =  str_extract(Sample_Name, "(?<=-)[^-]+$")) %>%
  merge(pool_size, by = 'species')
         
   
samples.as %>% head(3)



# Loading VCF ---------------------------------------------------------------------

# depth <- 2000
# Generate pooldata file
msas.raw <- vcf2pooldata(vcf.file = sprintf("../data/VCF/MSAS/SA/%s_raw_calls.d2000.filtered_calls.vcf", species),
                           poolsizes = samples.as[,"pool_size"], poolnames = samples.as[,"Sample_Name"])

msas.raw@poolnames


# Diagnostics
# Cases of zero-coverage
which(colSums(msas.raw@readcoverage) == 0)# %>% View()
# Calculate allele frequencies
frq.raw <- msas.raw@refallele.readcount/msas.raw@readcoverage 
colnames(frq.raw) <- msas.raw@poolnames


# Missing values

# Samples with high  number of NA
miss_count <- sapply(frq.raw %>% as.data.frame(), function(x) sum(is.na(x))) 
high_miss <- miss_count[which(miss_count > 1)] %>%as.data.frame() %>% rownames()
length(high_miss)
# Missing value rate 
length(which(is.na(frq.raw)))/length(frq.raw)



# return SNPs for min.cov.per.pool = 5
msas.lm <- pooldata.subset(msas.raw,  min.cov.per.pool = 5, verbose = FALSE, return.snp.idx = T)
# Calculate allele frequencies (high missing discarded)
frq <- msas.lm@refallele.readcount/msas.lm@readcoverage 
colnames(frq) <- msas.lm@poolnames


length(which(is.na(frq)))/length(frq)


# Compute PCA
#frq[is.na(frq)] <- 0
msas.pca <- prcomp(frq %>% t(), scale = T)




# Fst -----------------------------------------------------------------
# Genome-wide Fst all populations
msas.fst <- computeFST(msas.lm, nsnp.per.bjack.block = 5, method = "Identity")

# Pairwise Fst ------------------------------------------------------------



msas.pairwisefst <- compute.pairwiseFST(msas.lm, nsnp.per.bjack.block = 5, method = "Identity")
fstplot[[species]] <- plot_fstats(msas.pairwisefst)



msas.pairwisefst@values %>% head(3)



vis_prep <- function(df){
  
  
  df2 <- df %>% as.data.frame()
  colnames(df2) <- c("Fst", "Fst.mean", "se", "Q2", "Q2.mean", "Q2.se", "Nsnp")
  df2$comparison <- rownames(df)
  df2 <- df2 %>% 
    separate(comparison, sep = ";", into = c("S1", "S2")) %>%
    separate(S1, sep = "-", into = c("species1", "cultivar1", "rep1")) %>%
 
    separate(S2, sep = "-", into = c("species2", "cultivar2", "rep2")) %>%
  
  return(df2)
  
} 



cultivar_colors <- c("black", "orange3")#, "orange4", "brown")



msas.pairwisefst.viz1 <- vis_prep(msas.pairwisefst@values)
msas.pairwisefst.viz1 %>% head(3)

std <- function(x) sd(x)/sqrt(length(x))

r <- msas.pairwisefst.viz1 %>% mutate(comparison = ifelse(cultivar1 == cultivar2, 'Within','Between')) %>%
  group_by(comparison) %>% summarise(mean.Fst = mean(Fst), sd.Fst = sd(Fst), 
                                    # lower95 =  mean.Fst + c(-1.96)*se.Fst,
                                    # high95 =  mean.Fst + c(1.96)*se.Fst,
                                     fst.sd = sprintf("%.3f±%.3f", mean.Fst, sd.Fst)) %>% t()

results$mean.fst.between <- r[4] 
results$mean.fst.within <- r[8]

#png(filename = sprintf("05_Data_analysis/Figures/%s_within_between_Fst.png", species), height = 5, width = 3.75, res = 300, units = 'in')

plot.df <- msas.pairwisefst.viz1 %>% mutate(within = ifelse(cultivar1 == cultivar2, T,F)) %>%
  group_by(within) %>%
  mutate(within = ifelse(within, paste0("W (N=", n(), ")"), paste0("B (N=", n(), ")")),
         Fst    = ifelse(Fst < 0, 0, Fst)) %>%
  mutate( within = factor(within, levels=c(paste0("W (N=", n(), ")"), paste0("B (N=", n(), ")")) ) )
  # ggpubr::ggdotplot(x = "within", y = "Fst",   width = 0.45, size = 1, color = 'black', fill = NA, add.params = list(size = 4, alpha = 0.75)) + 


  
  
p <-   ggplot(data = plot.df, aes(x = within, y = Fst)) + 
  geom_jitter(fill = NA, size = 5, width = 0.2, shape = 21, alpha = 0.25) + 
  geom_boxplot(size = 0.5,  color = 'red', outlier.shape = NA) +
    
    
    # Add significance comparison
    stat_compare_means(aes(group = within), method = "wilcox", label = "p.signif") + 
    
    # Add horizontal grid lines for specific Fst values
    geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", color = "gray") + 
    
    # Title and labels
    ggtitle(species.df[species.df$species == species,]$fspecies) + 
    xlab("") + 
    
    # Customize theme
    theme(
      axis.title.y = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      plot.title = ggtext::element_markdown(),
      panel.background = element_rect(fill = "white", color = "white"),  # White background
      panel.grid.major = element_blank(), #(color = "gray", size = 0.2),  # Gray grid lines
      panel.grid.minor = element_blank(), # element_line(color = "gray", size = 0.2),
      legend.position = "none"  # Remove legend if not needed
    ) + 
    ylim(0, 1) + 
    coord_flip()
  
  
bxplots[[species]] <- p 



results$snp <- msas.lm@nsnp

results.final <- rbind(results.final, results)
}

# results.final %>% write.csv("../results/tables/fivespp_pure_Fst_stats.txt")


# Table 4 [SA samples] -------------------------------------------------------------

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

results.final %>% arrange(species) %>% # round_df(3) %>% 
  #mutate(mean.fst.95CI = sprintf("%.3f [%.3f, %.3f]", mean.fst,mean.fst.lower95CI,mean.fst.high95CI )) %>%
  #mutate(mean.fst.95CI.between = sprintf("%.3f [%.3f, %.3f]", mean.fst,mean.fst.lower95CI,mean.fst.high95CI )) %>%
  select(species, snp, mean.fst.between, mean.fst.within) %>% 
#  xtable::xtable() %>%
  write.table("../results/tables/MSAS/Table_4.SA_samples.txt", quote = F, sep = "\t")



# Pairwise Fst visualization ------------------------------------------------------------

msas.gbs.colors <- c("#E69F00", "#56B4E9", "#D55E00")

load("../data/topweide_results/pairwise_Fst.RData")
load("../data/topweide_results/GBS_pairwise_Fst.RData")

lp.ext.df <- lp.ext %>% 
  mutate(within = ifelse(within == 'Within\n(N=16)', 'W (N=16)','B (N=120)'),
         Fst    = ifelse(Fst < 0, 0, Fst)) %>% 
  mutate( within = factor(within, levels=c( "B (N=120)", "W (N=16)") ) )

levels(lp.ext.df$within)

lp.ext.p <-  ggplot(data = lp.ext.df, aes(x = within, y = Fst)) + 
  
  geom_jitter(fill = NA, size = 5, width = 0.2, shape = 21, alpha = 0.25) + 
  geom_boxplot(size = 0.5,  color = 'red' , outlier.shape = NA) +
  # Add significance comparison
  stat_compare_means(aes(group = within), method = "wilcox", label = "p.signif") + 
  
  # Add horizontal grid lines for specific Fst values
  geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", color = "gray") + 
  
  # Title and labels
  ggtitle("Extended *L. perenne* (MSAS)") + 
  xlab("") + 
  
  # Customize theme
  theme(
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = ggtext::element_markdown(),
    panel.background = element_rect(fill = "white", color = "white"),  # White background
    panel.grid.major = element_blank(), #(color = "gray", size = 0.2),  # Gray grid lines
    panel.grid.minor = element_blank(), # element_line(color = "gray", size = 0.2),
    legend.position = "none"  # Remove legend if not needed
  ) + 
  ylim(0, 1) + 
  coord_flip()
  
#lp.ext.p
  

  
  gbs.lp.ext.df <- gbs.lp.ext %>% 
  mutate(within = ifelse(within == 'Within\n(N=18)', 'W (N=18)','B (N=135)'),
         Fst    = ifelse(Fst < 0, 0, Fst)) %>%
  mutate( within = factor(within, levels=c(  "B (N=135)", "W (N=18)") ) )# %>%
  

  gbs.lp.ext.p <- ggplot(data = gbs.lp.ext.df, aes(x = within, y = Fst)) + 
    
    geom_jitter(fill = NA, size = 5, width = 0.2, shape = 21, alpha = 0.25) + 
    geom_boxplot(size = 0.5,  color = 'red', outlier.shape = NA) +
  # Add significance comparison
  stat_compare_means(aes(group = within), method = "wilcox", label = "p.signif") + 
  
  # Add horizontal grid lines for specific Fst values
  geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", color = "gray") + 
  
  # Title and labels
  ggtitle("Extended *L. perenne* (GBS)") + 
  xlab("") + 
  
  # Customize theme
  theme(
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = ggtext::element_markdown(),
    panel.background = element_rect(fill = "white", color = "white"),  # White background
    panel.grid.major = element_blank(), #(color = "gray", size = 0.2),  # Gray grid lines
    panel.grid.minor = element_blank(), # element_line(color = "gray", size = 0.2),
    legend.position = "none"  # Remove legend if not needed
  ) + 
  ylim(0, 1) + 
  coord_flip()
  
#  gbs.lp.ext.p



wilco.pairFst <- function(p){
r <- p + ylim(c(-0.12,1))  + 
  
  stat_compare_means(label = "p.signif", ref.group = "W\n(N=6)", method = "wilcox", size = 6, color = '#D55E00', label.y.npc = 0.8) +  
  theme(legend.position="none", panel.background = element_rect(fill = "gray95",
                                                                  colour = "gray95"))
return(r)
}

wilco.pairFst.lp <- function(p){
  r <- p + ylim(c(-0.12,1))  + 
    
    stat_compare_means(label = "p.signif", ref.group = "W\n(N=16)", method = "wilcox", size = 6, color = '#D55E00', label.y.npc = 0.8) +  
    theme(legend.position="none", panel.background = element_rect(fill = "gray95",
                                                                  colour = "gray95"))
  return(r)
}

wilco.pairFst.lp_gbs <- function(p){
  r <- p + ylim(c(-0.12,1))  + 
    
    stat_compare_means(label = "p.signif", ref.group = "W\n(N=18)", method = "wilcox", size = 6, color = '#D55E00', label.y.npc = 0.8) +  
    theme(legend.position="none", panel.background = element_rect(fill = "gray95",
                                                                  colour = "gray95"))
  return(r)
}

plot.clean_axis <- function(p){
  r <- p +  theme(axis.line.x=element_blank(), # axis.text.x.left =element_blank(),
                                           axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                                           axis.title.x=element_blank())
  return(r) }



bxplots2            <-  lapply( bxplots, plot.clean_axis)
bxplots2$Lp.ext     <- lp.ext.p %>% plot.clean_axis
bxplots2$Lp.ext_gbs <- gbs.lp.ext.p


p.panels <- cowplot::plot_grid(plotlist = bxplots2, nrow = 7, ncol = 1)


# saving Fig. 1 panel (a)-----------------------------------------------------------


ggsave(p.panels, filename="../results/figures/Fig_1.panel_A.raw.pdf", device="pdf", dpi = 300, width = 5, height = 8, units = 'in')

# saving Table 4 Extended L. perenne samples ----------------------------------

rbind( lp.ext %>% group_by(within, Nsnp) %>% summarise(mean.Fst = mean(Fst), sd.Fst = sd(Fst), 
                                              fst.sd = sprintf("%.3f±%.3f", mean.Fst, sd.Fst)) %>% #  %>% t()
            mutate(Dataset = "Extended L. perenne", Method = "MSAS"),
       
gbs.lp.ext %>% group_by(within, Nsnp) %>% summarise(mean.Fst = mean(Fst), sd.Fst = sd(Fst),
                                             fst.sd = sprintf("%.3f±%.3f", mean.Fst, sd.Fst))  %>% #  %>% t()
  mutate(Dataset = "Extended L. perenne", Method = "GBS")
) %>%
  write.table("../results/tables/extended_loper/Table_4.Ext_lper.txt", quote = F, sep = "\t")
