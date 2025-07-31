# MSAS Poolfstat analysis with pure and sp-sorted reads
# Author: Miguel Loera
# 02.04.2023

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
library(adegenet)
library(ggrepel)


# poolfstat ---------------------------------------------------------------

# species <- 'Tp'
results.final <- data.frame()
bxplots <- list()
permanova.list <- list()

species.df <- data.frame(species = c(
  #'Dg', 
  'Fp', 'Lp', 'Tp', 'Tr'),
  fspecies = c(
    # '*Dactylis<br>glomerata*<br>MSAS', 
    'T.A.<br>*Festuca<br>pratensis*<br>MSAS', 
    'T.A.<br>*Lolium<br>perenne*<br>MSAS', 
    'T.A.<br>*Trifolium<br>pratense*<br>MSAS', 
    'T.A.<br>*Trifolium<br>repens*<br>MSAS'))



# Define data path
#path_daten <- path_daten <- paste0(getwd(), "/") 
# Load sample names

species <- 'Tp'

pool_size = data.frame(species =  c('Dg', 'Fp', 'Lp', 'Tp', 'Tr'), pool_size = c(80,40,40,40,80))



for (species in species.df$species) {
  
  
  
  
  
  results <- data.frame(species = c(species))
  
  # sample names ------------------------------------------------------------
  
  
  samples <- read.csv(sprintf("../data/VCF/MSAS/SAMS/samples/%s_samples.txt", species), header = F, col.names = 'sample') %>%
    mutate(sample = gsub(".*/", "", sample),
           sample = gsub("\\..*", "", sample),
           sample = gsub("A-", "A100-", sample),
           sample = gsub("B-", "B100-", sample),
           sample2 = sample) %>%
    separate(sample2, sep = '-', into = c('mx', 'cultivar', 'rep'))
  
  
  
  
  # Loading mixed VCF ---------------------------------------------------------------------
  
  # depth <- 2000
  # Generate pooldata file
  msas.raw <- vcf2pooldata(vcf.file = sprintf("../data/VCF/MSAS/SAMS/%s_taxsortSNP.vcf", species),
                           poolsizes = rep(pool_size[which(species == species),]$pool_size,15), poolnames = samples[,"sample"])
  
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
  rownames(frq) <- paste(msas.lm@snp.info$Chromosome, 
                         msas.lm@snp.info$Position,
                         sep = "_")
  
  
  length(which(is.na(frq)))/length(frq)
  
  
  # Compute PCA
  #frq[is.na(frq)] <- 0
  msas.pca <- prcomp(frq %>% t(), scale = T)
  
  
  # Fst -----------------------------------------------------------------
  # Genome-wide Fst all populations
  msas.fst <- computeFST(msas.lm, nsnp.per.bjack.block = 5, method = "Identity")
  msas.fst$mean.fst
  
  results$snp <- msas.lm@nsnp
  # results$FST <- msas.fst$FST
  results$Fst <- msas.fst$Fst[1] %>% as.numeric() %>% round(3)
  # Pairwise Fst ------------------------------------------------------------
  
  msas.pairwisefst <- compute.pairwiseFST(msas.lm, nsnp.per.bjack.block = 5, method = "Identity")
  
  
  
  
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
  
  
  
 
  
  msas.pairwisefst.viz1 <- vis_prep(msas.pairwisefst@values)
  msas.pairwisefst.viz1 %>% head(3)
  
  std <- function(x) sd(x)/sqrt(length(x))
  
  r <- msas.pairwisefst.viz1 %>% filter(cultivar1 != 'AB50') %>% filter( cultivar2 != 'AB50') %>% mutate(comparison = ifelse(cultivar1 == cultivar2, 'Within','Between')) %>%
    group_by(comparison) %>% summarise(mean.Fst = mean(Fst), sd.Fst = sd(Fst), 
                                       # lower95 =  mean.Fst + c(-1.96)*se.Fst,
                                       # high95 =  mean.Fst + c(1.96)*se.Fst,
                                       fst.sd = sprintf("%.3f±%.3f", mean.Fst, sd.Fst)) %>% t()
  
  results$mean.fst.between <- r[4] 
  results$mean.fst.within <- r[8]
  

  
  results.final <- rbind(results.final, results)
  

# permanova ---------------------------------------------------------------

  perms = 9999
  frq.2 <- frq[,!colnames(frq) %in% c("Mx-AB50-1", "Mx-AB50-2", "Mx-AB50-3")]
  pops <- data.frame(samples = colnames(frq.2))
  
  pops <- pops %>% separate(samples, sep = '-', into = c('sp', 'pops', 'rep'))
  
  print(pops)
  
  permanova.rs <- vegan::adonis2(t(frq.2) ~ pops, data = pops, perm = perms, method = 'euclidean')
  print(permanova.rs)
  permanova.list[[species]] <- permanova.rs
  hist(frq, title = species)  
  
}




# Table 4 SAMS samples ----------------------------------------------------



results.final |> write.table("../results/tables/MSAS/Table_4.SAMS_samples.txt", sep = "\t", quote = F, row.names = F)




# Table S3 PERMANOVA SAMS samples ---------------------------------------------------


nombres <- names(permanova.list)
do.call(rbind, lapply(nombres |> seq_along(), function(i){nombre = nombres[i]; print(nombre); df <- permanova.list[[nombre]] |> as.data.frame() |> mutate(sp = nombre |> as.character() ); return(df) }))
