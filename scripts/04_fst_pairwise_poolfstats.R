
# loading libraries -------------------------------------------------------

library(tidyverse)


# loading data ------------------------------------------------------------

gbs.fst_raw           <- read.csv('../results/tables/extended_loper/fst_pairwise_gbs.csv', row.names = 1)
colnames(gbs.fst_raw) <- gsub(colnames( gbs.fst_raw ), pattern = '\\.', replacement = '-')

gbs.fst_raw %>% head(3)


msas.fst_raw           <- read.csv('../results/tables/extended_loper/fst_pairwise_msas.csv', row.names = 1)
colnames(msas.fst_raw) <- gsub(colnames( msas.fst_raw ), pattern = '\\.', replacement = '-')

msas.fst_raw %>% head(3)


# Visualization prep. -------------------------------------------------------------

gbs.fst_raw.2 <- gbs.fst_raw %>% dplyr::mutate(s1 = rownames(gbs.fst_raw)) %>%
  reshape2::melt(id.vars = c('s1')) %>%
  dplyr::rename(s2 = variable, Fst.gbs = value) %>%
  dplyr::mutate(c1 = substr(s1, 1, 3), c2 = substr(s2, 1, 3)
                )


msas.fst_raw.2 <- msas.fst_raw %>% dplyr::mutate(s1 = rownames(msas.fst_raw)) %>%
  reshape2::melt(id.vars = c('s1')) %>%
  dplyr::rename(s2 = variable, Fst.msas = value) %>%
  dplyr::mutate(c1 = substr(s1, 1, 3), c2 = substr(s2, 1, 3)
  )

gbs_msas.fst <- merge(gbs.fst_raw.2, msas.fst_raw.2, by = c('s1', 's2', 'c1', 'c2')) %>%
  dplyr::filter( s1 != s2 ) %>%
  dplyr::mutate(

    comparison_type = ifelse(
      c1 == c2, 'Within', 'Between' ))


# Fig. 2 Panel (c) -----------------------------------------------------------

gbs_msas.fst.2  <- gbs_msas.fst %>%
  dplyr::filter( s1 != s2 ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    
    comparison = paste(  sort( c(s1 %>% as.character() , s2 %>% as.character() ) )[1],  sort( c(s1 %>% as.character() , s2 %>% as.character() ) )[2] )
    
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(comparison) %>%
  dplyr::mutate(
    comparison_num = seq_along(comparison)  # Add a consecutive number for each comparison
  ) %>%
  dplyr::filter( comparison_num == 1) 


# Fit the linear model
lm_model <- lm(Fst.gbs ~ Fst.msas, data = gbs_msas.fst.2)

# Calculate residuals
gbs_msas.fst.2$fitted_values <- predict(lm_model)
gbs_msas.fst.2$residuals <- residuals(lm_model)


# Get the model summary
model_summary <- summary(lm_model)

# Extract coefficients, R-squared, and p-value
intercept <- coef(model_summary)[1]
slope <- coef(model_summary)[2]
r_squared <- model_summary$r.squared
p_value <- model_summary$coefficients[2, 4]  # Extracting p-value from the coefficients matrix
# Create a formula string for annotation
formula_label <- paste0("y = ", round(intercept, 3), " + ", round(slope, 3), " * x")
stats_label <- paste0("R² = ", round(r_squared, 3), ", p = ", format.pval(p_value, digits = 2, format = "f"))

Fst.limit <- gbs_msas.fst.2 %>% dplyr::filter( comparison_type == 'Within') %>% dplyr::ungroup() %>% dplyr::reframe(max.Fst = max(Fst.gbs, na.rm = T)) %>% as.numeric()

plot5.2 <- ggplot(gbs_msas.fst.2) +
  annotate("text", x = Inf, y = Inf, label = formula_label, hjust = 1.1, vjust = 2, size = 5, color = "red") +
  annotate("text", x = Inf, y = Inf, label = stats_label, hjust = 1.1, vjust = 3.5, size = 5, color = "red") +
  # # Plot residual lines
  # geom_segment(aes(x = Fst.msas, y = Fst.gbs, xend = Fst.msas, yend = fitted_values), 
  #              color = "gray80", size = 1) +  # Residuals as blue lines
  # add hline at gbs.fst = 0.04
  geom_hline( yintercept = Fst.limit, linetype = 2, size = 1, color = 'black') +
  # Plot original data points
  geom_point(aes(x = Fst.msas, y = Fst.gbs, col = comparison_type), size = 6, shape = 21) +
  scale_color_manual( values = c('#E69F00', '#56B4E9')) +
  # Plot regression line
  geom_smooth(aes(x = Fst.msas, y = Fst.gbs), method = "lm", se = FALSE, color = "red", size = 2, inherit.aes = F, alpha = 0.2) +
  

  
  labs(# title = "Regression Plot with Residual Lines",
       x = "Pairwise Fst (MSAS)",
       y = "Pairwise Fst (GBS)") +
  
  #theme_minimal() +
  theme( legend.position = 'top')

plot5.2

ggsave(plot = plot5.2, filename  = "../results/figures/Fig_2.panel_C.raw.pdf", width = 8, height = 8, dpi = 300, units = 'in')

