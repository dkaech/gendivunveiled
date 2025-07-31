# Taxonomic assignment rate
# (calculated based on mapping rates to msas_ref.fasta)

library(tidyverse)
library(ggplot2)
library(reshape2)

# loading data ------------------------------------------------------------


data.df <- read.csv("../data/VCF/MSAS/taxsort_idxstats_f2_F256.txt", sep = '\t', header = F, 
                    col.names = c('SAMPLE', 'amplicon', 'length', 'reads', 'unmapped')) %>%
  separate(amplicon, sep = '_', into = c('ASSP', 'AMPLICON'))   

data.df$SAMPLE <- gsub(data.df$SAMPLE, pattern = 'TAXSORT_BAM/', replacement = '')
data.df$SAMPLE <- gsub(data.df$SAMPLE, pattern = '.msas_ref.sam', replacement = '')

data.df2 <- data.df %>%
  mutate(REFSP = str_sub(SAMPLE, 1, 2), REFCULT = str_sub(SAMPLE, 1,4)) %>%
  filter(REFSP != 'Mx' & ASSP != '*')

# summaries of mapped and unmapped reads ---------------------------------------------------------------

data.df$reads %>% sum()
data.df$unmapped %>% sum()
# summaries of tax. assignment --------------------------------------------------------------------
# tax. assignment by sample in percent reads
t1 <- data.df2 %>%
  group_by(SAMPLE, REFSP, ASSP) %>%
  summarise(reads = sum(reads)) %>%
  group_by(SAMPLE) %>% 
  mutate(scaled_reads = reads/sum(reads), sum_reads = sum(reads)) %>% #%>% # View()
  select(SAMPLE, REFSP, ASSP, reads, scaled_reads )

# tax. assignment by species in percent reads
t1.2 <- data.df2 %>%
  group_by(REFSP, ASSP) %>%
  summarise(reads = sum(reads)) %>%
  group_by(REFSP) %>% 
  mutate(scaled_reads = reads/sum(reads), sum_reads = sum(reads)) %>% #%>% # View()
  select(REFSP, ASSP, reads, scaled_reads )

  # tax. assignment at > 1% of reads per species
  t1.2 %>% filter(scaled_reads > 0.01)
  


# visualization: assignment by species in percent reads -----------------------------------------------------------

  
  # boxplot of tax. assignment by species in percent reads
  
  p1 <- t1 %>% #filter(scaled_reads > 0) %>%
   ggplot(aes(x = ASSP, y = scaled_reads)) + 
    facet_wrap(REFSP ~ .) +
    geom_boxplot() +
    geom_jitter(alpha = 0.3, col = 'blue', width = 0.1) +
    xlab('Assigned species') +
    ylab('# of reads') +
    ggthemes::theme_clean()
  
  p1

# visualizations: taxonomical assignment rate ----------------------------------------------------------
  
  # boxplot of read assignment by species
  p2 <- data.df2 %>%
    group_by(SAMPLE, REFSP) %>%
    mutate( total_reads = sum(reads), match = ifelse(REFSP == ASSP, 'Match', 'Missmatch')) %>%
    group_by(SAMPLE, REFSP, match) %>%
    summarise(rate = sum(reads)*100/unique(total_reads) ) %>%
    ggplot(aes(x = match, y = rate)) + 
    facet_wrap(REFSP ~ ., ncol= 5) +
    geom_boxplot() +
    geom_jitter(alpha = 0.3, col = 'blue', width = 0.1) +
    ggthemes::theme_clean() +
    xlab("") + ylab("Assignment rate") + theme(legend.position="none", panel.background = element_rect(fill = "gray95",
                                                                                                       colour = "gray95"))
  # table of read assignment rate by species
  data.df2 %>%
    group_by(REFSP) %>%
    mutate( total_reads = sum(reads), match = ifelse(REFSP == ASSP, 'Match', 'Missmatch')) %>%
    group_by(REFSP, match) %>%
    reframe(rate = sum(reads)*100/unique(total_reads) ) %>%
    write.table("../results/tables/MSAS/additional_table.Taxonomic_assignment_rate.MS_samples.txt", sep = "\t", quote = F, row.names = F)
  
  # boxplot of read assignment rate by accession
  
  data.df2 %>%
    group_by(SAMPLE, REFCULT) %>%
    mutate( total_reads = sum(reads), match = ifelse(REFSP == ASSP, 'Match', 'Missmatch')) %>%
    group_by(SAMPLE, REFCULT, match) %>%
    summarize(rate = sum(reads)*100/unique(total_reads) ) %>%
    ggplot(aes(x = match, y = rate)) + 
    facet_wrap(REFCULT ~ .) +
    geom_boxplot()


# visualizations: taxonomical assignment rate & summary read stats by species --------------------------------------------------------------------


data.df3 <- data.df %>%
  mutate(REFSP = str_sub(SAMPLE, 1, 2), REFCULT = str_sub(SAMPLE, 1,4)) %>%
  filter(REFSP != 'Mx') %>%
  
  group_by(REFSP) %>%
  mutate( total_reads = sum(reads), match = ifelse(REFSP == ASSP, 'match', 'missmatch')) %>%
  group_by(REFSP, match) %>%
  summarize(rate = sum(reads)*100/unique(total_reads) )  %>%
  dcast(REFSP ~ match, value.var = 'rate', fun.aggregate = sum)

  
data.df %>%
  mutate(REFSP = str_sub(SAMPLE, 1, 2), REFCULT = str_sub(SAMPLE, 1,4)) %>%
  filter(REFSP != 'Mx') %>%
  group_by(REFSP) %>%
  mutate( total_reads = sum(reads)) %>%
  group_by(REFSP, ASSP) %>%
  summarize(rate = sum(reads)*100/unique(total_reads), n.reads = unique(total_reads) ) %>%
  merge(data.df3, by = 'REFSP' ) %>%
  ggplot(aes(x = ASSP, y = rate)) +
  geom_text(aes(label = sprintf('Total reads: %d\nmatch: %.2f%%\nmissmatch: %.2f%%', n.reads, match, missmatch)), x = 1, y = 80, hjust = 0, inherit.aes = F, color = 'red' ) +
  facet_grid(REFSP ~ .) +
  geom_bar(stat = 'identity') +
  ggthemes::theme_clean() +
  xlab('Assigned species') +
  ylab('Assignment rate [%]') +
  theme( panel.spacing.x = unit(1,"line"),  
         panel.spacing.y = unit(1.5,"line"))
