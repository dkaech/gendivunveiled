# Purpose: Load required packages and define paths
# Author: Damian Käch <dkaech@ethz.ch>

# delete all objects from environment
rm(list = ls())

cat("\n\nAll objects deleted!\n\n")
cat("--------------------\n\n")


# Path definition ---------------------------------------------------------

# Path VCF files
path_vcf <- "../data/VCF/"
# Path figures 
path_figures <- "../results/figures/extended_loper/"
# Path tables
path_tables <- "../results/tables/extended_loper/"
# Path admixture GBS
path_admix <- "../data/ADMIXTURE/GBS/"
# Path structure GBS
path_strct <- "../data/STRUCTURE/GBS/"
# Path admixture MSAS
path_admix_ampse <- "../data/ADMIXTURE/MSAS/"
# Path structure MSAS
path_strct_ampse <- "../data/STRUCTURE/MSAS/"

# Package loading ---------------------------------------------------------

# Load packages
cat("Loading packages...\n")
# Required packages
pkgs <- c("dplyr", "tidyr", "stringr", "ggplot2", "poolfstat", "adegenet", "ggpubr", "ggrepel", "purrr")

# Get the names installed packages
inst_pkgs <- installed.packages()[, "Package"]

# Install required packages from CRAN
req_pkgs <- setdiff(pkgs, inst_pkgs)
if (length(req_pkgs) > 0) {
  install.packages(req_pkgs, clean = TRUE)
}

# Load necessary packages
suppressMessages(lapply(c(pkgs), library, character.only = TRUE))

# Save current session info
session_info <- sessionInfo()

# Remove obsolete variables
rm(inst_pkgs, req_pkgs)

cat("\nPackages loaded!\n\n")
cat("--------------------\n\n")


# Function definition -----------------------------------------------------

# Get upper and lower triangle of correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Pairplot
pairplot <- function(dat_pairplot, poly_order = 1, custom_alpha = 1, 
                     save = FALSE, file_type = "pdf", 
                     plot_height = 200, plot_width = 350, plot_resolution = 400) {
  # Define new lists for each plot type 
  list_gg_hist <- list()
  list_gg_point <- list()
  list_gg_numb <- list()
  
  # Generate plots and write into list
  for (colnumber in 1:ncol(dat_pairplot)) {
    
    colname <- colnames(dat_pairplot)[colnumber]
    
    # Histogram
    gg_hist <- ggplot(dat_pairplot, aes_string(colname)) +
      geom_histogram(fill = "#D55E00", color = "black") +
      ggtitle(colname) +
      theme_bw() +
      theme(plot.title = element_text(color = "black", hjust = 0.5, size = 30),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank())
    
    list_gg_hist[[colname]] <- gg_hist
    
    for (colnumber_2 in (colnumber+1):ncol(dat_pairplot)) {
      
      if (colnumber_2 <= ncol(dat_pairplot)) {
        
        colname_2 <- colnames(dat_pairplot)[colnumber_2]
        
        # Correlation number
        cor_number <- round(cor(dat_pairplot[,colname],dat_pairplot[,colname_2], 
                                use = "complete.obs", method = "spearman"),
                            digits = 2)
        # text_size <- 5 + cor_number^2 * 150/ncol(dat_pairplot)
        text_size <- cor_number * plot_height / ncol(dat_pairplot)
        gg_numb <- data.frame(x = 0.5, y = 0.5, text = as.numeric(cor_number)) %>% 
          ggplot() +
          # Text size depending on cor_number, ^2 -> positive value
          geom_text(aes(x, y, label = text), size = text_size) +
          theme_bw() +
          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank())
        
        list_gg_numb[[paste(colname,":",colname_2)]] <- gg_numb
      }
    }
    
    for (colnumber_3 in 1:colnumber) {
      if (colnumber != colnumber_3) {
        colname_3 <- colnames(dat_pairplot)[colnumber_3]
        # Scatterplot
        gg_point <- ggplot(dat_pairplot, aes_string(colname_3, colname)) +
          geom_point(shape = 1, alpha = custom_alpha) +
          geom_smooth(method = "lm", color = "#D55E00", formula = y ~ poly(x, poly_order)) +
          theme_bw() +
          theme(##axis.text = element_blank(),
            axis.title = element_blank(),
            #axis.ticks = element_blank()
          )
        
        list_gg_point[[paste(colname_3,":",colname)]] <- gg_point
      }
    }
  }
  # Delete last element of the lists as it is the interaction with itself
  list_gg_numb <- list_gg_numb[-length(list_gg_numb)]
  
  # Transform the three lists into one large list so that the 
  # order is suitable for ggarrange (fills plots per row from left to right)
  list_pair_plot <- list()
  for (n in 0:(ncol(dat_pairplot)-1)) {
    point <- n
    hist <- n+1
    numb <- ncol(dat_pairplot)-1-n
    
    element_point <- list_gg_point[0:point]
    element_hist <- list_gg_hist[hist]
    element_numb <- list_gg_numb[0:numb]
    
    if (length(element_point)>0){
      listlength <- length(list_pair_plot)
      list_pair_plot[c((listlength + 1):(listlength + length(element_point)))] <- element_point
    }
    
    listlength <- length(list_pair_plot) + 1
    list_pair_plot[listlength] <- element_hist
    
    if (length(element_numb)>0){
      listlength <- length(list_pair_plot)
      list_pair_plot[c((listlength + 1):(listlength + length(element_numb)))] <- element_numb
    }
    
    # Delete used elements
    if (length(element_point)>0) {
      list_gg_point <- list_gg_point[-c(1:length(element_point))]
    }
    if (length(element_numb)>0) {
      list_gg_numb <- list_gg_numb[-c(1:length(element_numb))] 
    }
  }
  cat("-----------\n\n")
  cat("Plotting...\n\n")
  gg_pairplot <- ggarrange(plotlist = list_pair_plot,
                           ncol = ncol(dat_pairplot),
                           nrow = ncol(dat_pairplot)) %>% suppressMessages()
  print(gg_pairplot)
  if (save) {
    # Ask user for file name
    question <- NULL
    repeat {
      filename <- readline(prompt = message("\n-----------\n\nType in the file name under which you want to save the plot:\n"))
      if (str_count(filename) > 0) {
        break
      }
    } 
    # Save plot under the user's input
    cat("Saving...\n\n")

    ggsave(filename = paste0(getwd(), "/", filename, ".", file_type), 
           plot = last_plot(), units = "mm", #limitsize = FALSE,
           height = plot_height, width = plot_width, 
           dpi = plot_resolution, device = file_type)
    
    cat("Plot saved in your working directory under:\n\n", paste0(filename, ".", file_type, "\n\n"))
    cat("-----------------------------------\n\n")
    
  }
}
