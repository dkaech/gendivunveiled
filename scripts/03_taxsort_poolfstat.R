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
library(ggplotify)


# sample names ------------------------------------------------------------


samples <- read.csv("../data/VCF/MSAS/mx_samples.txt", header = F, col.names = 'sample2') %>%
  mutate(sample = sample2) %>%
  separate(sample2, sep = '-', into = c('mx', 'cultivar', 'rep'))


# poolfstat ---------------------------------------------------------------


results.final <- data.frame()
bxplots <- list()

species.df <- data.frame(species = c(
  #'Dg', 
                                     'Fp', 'Lp', 'Tp', 'Tr'),
                         fspecies = c(
                          # '*Dactylis<br>glomerata*<br>MSAS', 
                           '*Festuca<br>pratensis*<br>MSAS', 
                           '*Lolium<br>perenne*<br>MSAS', 
                           '*Trifolium<br>pratense*<br>MSAS', 
                           '*Trifolium<br>repens*<br>MSAS'))



pool_size = data.frame(species =  c('Dg', 'Fp', 'Lp', 'Tp', 'Tr'), pool_size = c(80,40,40,40,80))



for (species in species.df$species) {
  
  results <- data.frame(species = c(species))
  
  
  
  
  
  # Loading VCF ---------------------------------------------------------------------
  
  # depth <- 2000
  # Generate pooldata file
  msas.raw <- vcf2pooldata(vcf.file = sprintf("../data/VCF/MSAS/MS/%s_taxsorted.vcf", species),
                           poolsizes = rep(pool_size[which(species == species),]$pool_size,9), poolnames = samples[,"sample"])
  
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
  results$Fst <- msas.fst$Fst[1] %>% as.numeric() %>% round(3)
  # Pairwise Fst ------------------------------------------------------------
  
  
  
  msas.pairwisefst <- compute.pairwiseFST(msas.lm, nsnp.per.bjack.block = 5, method = "Identity")
  
  

  
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



  
  msas.pairwisefst.viz1 <- vis_prep(msas.pairwisefst@values)
 
  std <- function(x) sd(x)/sqrt(length(x))
  
  r <- msas.pairwisefst.viz1 %>% filter(cultivar1 != 'AB50') %>% filter( cultivar2 != 'AB50') %>% mutate(comparison = ifelse(cultivar1 == cultivar2, 'Within','Between')) %>%
    group_by(comparison) %>% summarise(mean.Fst = mean(Fst), sd.Fst = sd(Fst), 
                                       # lower95 =  mean.Fst + c(-1.96)*se.Fst,
                                       # high95 =  mean.Fst + c(1.96)*se.Fst,
                                       fst.sd = sprintf("%.3f±%.3f", mean.Fst, sd.Fst)) %>% t()
  
  results$mean.fst.between <- r[4] 
  results$mean.fst.within <- r[8]
  
  
  results$snp <- msas.lm@nsnp
  
  results.final <- rbind(results.final, results)
}


# Table 4 MS samples -------------------------------------------------------------



results.final %>% arrange(species) %>% # round_df(3) %>% 
  
  select(species, snp, Fst, mean.fst.between, mean.fst.within) %>% 
  write.table("../results/tables/MSAS/Table_4.MS_samples.txt", quote = F, sep = "\t", row.names = F)
  #xtable::xtable()
