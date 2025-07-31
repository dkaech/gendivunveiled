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
library(adegenet)
library(ggrepel)

# sample names ------------------------------------------------------------


samples <- read.csv("../data/VCF/MSAS/mx_samples.txt", header = F, col.names = 'sample2') %>%
  mutate(sample = sample2) %>%
  separate(sample2, sep = '-', into = c('mx', 'cultivar', 'rep'))

# poolfstat ---------------------------------------------------------------

# species <- 'Tp'
results.final <- data.frame()
bxplots <- list()
dapc.ggp <- list()
species.df <- data.frame(species = c(
  #'Dg', 
  'Fp', 'Lp', 'Tp', 'Tr'),
  fspecies = c(
    # '*Dactylis<br>glomerata*<br>MSAS', 
    '*Festuca<br>pratensis*<br>MSAS', 
    '*Lolium<br>perenne*<br>MSAS', 
    '*Trifolium<br>pratense*<br>MSAS', 
    '*Trifolium<br>repens*<br>MSAS'))


species.df <- data.frame(species = c(
  #'Dg', 
  'Fp', 'Lp', 'Tp', 'Tr'),
  fspecies = c(
    # '*Dactylis glomerata* MSAS', 
    '*Festuca pratensis* MSAS', 
    '*Lolium perenne* MSAS', 
    '*Trifolium pratense* MSAS', 
    '*Trifolium repens* MSAS'))

# species <- 'Dg'

pool_size = data.frame(species =  c('Dg', 'Fp', 'Lp', 'Tp', 'Tr'), pool_size = c(80,40,40,40,80))



for (species in species.df$species) {
  

  # Loading taxsorted MS samples VCF ---------------------------------------------------------------------
  

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
   
  
  # Allele frequencies: taxsorted MS samples ------------------------------------------
  
    
     # Calculate allele frequencies (high missing discarded) -> taxsorted samples
    frq <- msas.lm@refallele.readcount/msas.lm@readcoverage 
    colnames(frq) <- msas.lm@poolnames
    rownames(frq) <- paste(msas.lm@snp.info$Chromosome, 
                           msas.lm@snp.info$Position,
                               sep = "_")
    
    
    length(which(is.na(frq)))/length(frq)
  
  
  # Loading pure VCF (SA samples) ---------------------------------------------------------------------
  
  samples.as <- data.frame(Sample_Name = c(paste0(sprintf("%s-A-", species), 1:3), paste0(sprintf("%s-B-", species), 1:3)))
  
  samples.as <- samples.as %>% 
    mutate(species = str_extract(Sample_Name, "^[^-]+"),
           cultivar = str_extract(Sample_Name, "(?<=-)[^-]+"),
           rep =  str_extract(Sample_Name, "(?<=-)[^-]+$")) %>%
    merge(pool_size, by = 'species')
  
  
  samples.as %>% head(3)
  
  
  # depth <- 2000
  # Generate pooldata file
  msas_pure.raw <- vcf2pooldata(vcf.file = sprintf("../data/VCF/MSAS/SA/%s_raw_calls.d2000.filtered_calls.vcf", species),
                           poolsizes = samples.as[,"pool_size"], poolnames = samples.as[,"Sample_Name"])
  
  msas_pure.raw@poolnames
  
  
  # Diagnostics
  # Cases of zero-coverage
  which(colSums(msas_pure.raw@readcoverage) == 0)# %>% View()
  
  # Allele frequencies: SA samples ------------------------------------------
  
    
    # Calculate allele frequencies
    frq_pure.raw <- msas_pure.raw@refallele.readcount/msas_pure.raw@readcoverage 
    colnames(frq_pure.raw) <- msas_pure.raw@poolnames
    
    
    # Missing values
    
    # Samples with high  number of NA
    miss_pure_count <- sapply(frq_pure.raw %>% as.data.frame(), function(x) sum(is.na(x))) 
    high_pure_miss <- miss_pure_count[which(miss_pure_count > 1)] %>%as.data.frame() %>% rownames()
    length(high_pure_miss)
    # Missing value rate 
    length(which(is.na(frq_pure.raw)))/length(frq_pure.raw)
    
    
    
    # return SNPs for min.cov.per.pool = 5
    msas_pure.lm <- pooldata.subset(msas_pure.raw,  min.cov.per.pool = 5, verbose = FALSE, return.snp.idx = T)
    # Calculate allele frequencies (high missing discarded) -> SA samples
    frq_pure <- msas_pure.raw@refallele.readcount/msas_pure.raw@readcoverage 
    colnames(frq_pure) <- msas_pure.raw@poolnames
    rownames(frq_pure) <- paste(msas_pure.raw@snp.info$Chromosome, 
                                msas_pure.raw@snp.info$Position,
                           sep = "_")
    
    # subsetting SA sample allele frequencies for SNPs shared with MS samples
    frq_pure_ft <- subset(frq_pure, rownames(frq_pure) %in% rownames(frq))
   
  # DAPC --------------------------------------------------------------------
  # Taking the allele frequencies for shared SNPs between MS and SA samples
    
  dapc_all_mx <- bind_rows(
    t(frq_pure_ft) %>% as.data.frame(),  ## Allele frequencies in SA samples
    
    t(subset(frq, rownames(frq) %in% rownames(frq_pure_ft))  ## Allele frequencies in MS samples
      ) %>% 
      
      as.data.frame() ) %>% 
      mutate(Sample = rownames(.),
             Group = gsub("-.$", "", Sample) %>% 
               factor(),
             Group_nr = as.numeric(Group)) %>% 
      select(Sample, Group, Group_nr, contains("_"))
  
  snps <- colnames(dapc_all_mx %>% select(contains(paste0(species,"_"))))  
   
  dapc_all <- dapc(dapc_all_mx %>%
                     # Keep only numerical variables
                     select(-c("Sample","Group","Group_nr")),  grp = dapc_all_mx$Group, n.da = 5, n.pca = 4) 
  scatter(dapc_all, scree.pca = T)
  
    
    # DAPC plot ---------------------------------------------------------------
    df_grp2_1 <- dapc_all$grp.coord %>% 
      as.data.frame() %>% 
      mutate(Name = rownames(.),
             Category = "Group")
    
    df_grp2_2 <- dapc_all$ind.coord %>% 
      as.data.frame()  %>% 
  
      mutate(Category = "Ind") # %>% 
      
    
    df_grp2 <- bind_rows(df_grp2_1, df_grp2_2) %>% 
      mutate(Colors_3 = case_when(str_detect(rownames(.), paste0(species, "-A")) ~ "PA",
                                  str_detect(rownames(.), paste0(species, "-B")) ~ "PB",
                                  str_detect(rownames(.), "A100") ~ "A100",
                                  str_detect(rownames(.), "AB50") ~ "AB50",
                                  TRUE ~ "B100"))
    
    
    
    p2 <- ggplot(data = df_grp2, aes(x = LD1, y = LD2)) +
      geom_point(aes(shape = Category, 
                     color = Colors_3, 
                     alpha = Category,
                     size = Category)) +
      
      
      scale_shape_manual(values = c("Group" = 4,
                                    #"Individual_pure" = 19,
                                    "Ind" = 15)) +
      scale_alpha_manual(values = c("Group" = 1,
                                    #"Individual_pure" = 0.3,
                                    "Ind" = 0.3)) +
      scale_size_manual(values = c("Group" = 5,
                                   #"Individual_pure" = 10,
                                   "Ind" = 10)) +
   
      geom_text_repel(data = df_grp2, 
                      aes(label = Name, color = Colors_3),
                      max.overlaps = 50, 
                      # Set distance between labels and points
                      box.padding = 1,
                      # Bellow following segment length, no
                      # segments are drawn (we don't want segments)
                      min.segment.length = 2,
                      size = 5.3) +
      scale_color_manual(values = c("A100" = "blue4",
                                    "PA" = "#56B4E9",
                                    "PB" = "#D55E00",
                                    "B100" = "orange4",
                                    "AB50" = "green2")) +
      geom_vline(xintercept = 0, linetype = 3) +
      geom_hline(yintercept = 0, linetype = 3) +
      theme_bw() +
      theme(legend.position = "none", 
            axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            panel.background = element_rect(fill='transparent'), 
            plot.background = element_rect(fill='transparent', color = NA)
      )
    
    
    dapc.ggp[[species]] <- p2 + ggtitle(species.df[species.df$species == species,]$fspecies) +
      theme(axis.title.y=element_text(size=12), axis.text.y=element_text(size=12),
            axis.title.x=element_text(size=12), axis.text.x=element_text(size=12), plot.title = ggtext::element_markdown()) +
      ylim(-15,15) + xlim(-30,30)
    
    
    }
    

# Figure 1 panel (b) ------------------------------------------------------


  cowplot::plot_grid(plotlist = dapc.ggp) %>% ggsave(filename = "../results/figures/Fig_1.panel_B.raw.pdf", height = 12, width = 12, dpi = 'print')
