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



# Loading data & permanova ------------------------------------------------------------
set.seed(19900829)
results.final <- data.frame()
bxplots <- list()
permanova.list <- list()
species.df <- data.frame(species = c('Dg', 'Fp', 'Lp', 'Tp', 'Tr'),
                         fspecies = c(
                           '*Dactylis<br>glomerata*<br>MSAS', 
                           '*Festuca<br>pratensis*<br>MSAS', 
                           '*Lolium<br>perenne*<br>MSAS', 
                           '*Trifolium<br>pratense*<br>MSAS', 
                           '*Trifolium<br>repens*<br>MSAS'))



# Define data path
#path_daten <- path_daten <- paste0(getwd(), "/") 
# Load sample names

# species <- 'Dg'

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
msas.lm <- pooldata.subset(msas.raw, verbose = FALSE, return.snp.idx = T)
# Calculate allele frequencies (high missing discarded)
frq <- msas.lm@refallele.readcount/msas.lm@readcoverage 
colnames(frq) <- msas.lm@poolnames

# frq %>% View()
length(which(is.na(frq)))/length(frq)


# Compute PCA
#frq[is.na(frq)] <- 0
msas.pca <- prcomp(frq %>% t(), scale = T)



# permanova ---------------------------------------------------------------

perms = 9999

pops <- data.frame(samples = msas.lm@poolnames)

pops <- pops %>% separate(samples, sep = '-', into = c('sp', 'pops', 'rep'))

print(pops)

permanova.rs <- vegan::adonis2(t(frq) ~ pops, data = pops, perm = perms, method = 'euclidean')
#print(permanova.rs)
permanova.list[[species]] <- permanova.rs
hist(frq, title = species)


}


# Permanova results -------------------------------------------------------


nombres <- names(permanova.list)
do.call(rbind, lapply(nombres |> seq_along(), function(i){nombre = nombres[i]; print(nombre); df <- permanova.list[[nombre]] |> as.data.frame() |> mutate(sp = nombre |> as.character() ); return(df) }))
