#############################################################################
##
## Analysis script for the paper "Dynamic MR of muscle contraction during 
## electrical muscle stimulation as a potential diagnostic tool for 
## neuromuscular disease"
##
## Copyright 2024 Francesco Santini <francesco.santini@unibas.ch>
##
## Released under an Apache 2.0 license. See LICENSE for details
##
## Execution notes: does not work with Rscript. Execute under Rstudio.
## R version 4.1.2
##
############################################################################

dat <- read.csv('results_dyn_quant.csv')
force_dat <- read.csv('force.csv')

fig_dir <- 'figs/'
strain_to_use <- 'strain_1'

col2="skyblue2"
col1="mediumpurple3"
library(ggpubr)
library(ggplot2)
library(extrafont)
library(dplyr)
library(tidyr)
loadfonts()

###############################
##
## Data preparation
##
###############################

dat[dat$variable == "buildup","value"] <- - dat[dat$variable == "buildup","value"]
dat[dat$variable == "release","value"] <- - dat[dat$variable == "release","value"]

new_rows <- dat %>%
  filter(variable == "FF") %>%
  group_by(subject, pat_id, subj_type, pathology_type, variable) %>%
  summarize(value = mean(value), .groups = "drop") %>%
  mutate(ROI = "Total")

dat <- bind_rows(dat, new_rows)

new_rows <- dat %>%
  filter(variable == "t2") %>%
  group_by(subject, pat_id, subj_type, pathology_type, variable) %>%
  summarize(value = mean(value), .groups = "drop") %>%
  mutate(ROI = "Total")

dat <- bind_rows(dat, new_rows)

# Create new rows for FF100
new_rows_ff100 <- dat %>%
  filter(variable == "FF") %>%
  mutate(variable = "FF100", value = round(value * 100, digits=1))

# Bind the new rows to the original data frame
dat <- bind_rows(dat, new_rows_ff100)


# sane data, in wide format
df_wide <- dat %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  filter(subj_type %in% c("V", "P")) # Ensure we only have V and P types if necessary

# Binarize subject for logistic regression
df_wide$subj_binary <- as.integer(df_wide$subj_type == 'P')

###############################
##
## Helper functions
##
###############################

save_plot <- function(g) {
  ggsave(g, file= sprintf("%s%s%s",fig_dir, g$labels$title, ".svg"), width=6, height=6)
  ggsave(g, file= sprintf("%s%s%s",fig_dir, g$labels$title, ".png"), width=6, height=6)
  g
}

ggplot_theme <- theme(panel.background = element_rect(fill = 'white',colour='gray'))+theme(panel.grid.major = element_line(colour = "gray", size = 0.2), panel.grid.major.x = element_blank()) +
  theme(text=element_text(size=16, family="Arial"), plot.title = element_text(size=20, family="Arial"))

make_plot <- function(variable, ylab_str, title) {
  g <- ggplot (aes(x= subj_type, y=value, label=subject),data=dat[dat$variable==variable&dat$ROI=='Total',])+ylab(ylab_str)+xlab('Subjects')+
    geom_boxplot(alpha=c(0.2,0.2),col=c(col1,col2),fill=c(col1,col2), lwd=0.8, width=0.3)+
    geom_point(col="gray", size=1.5)+ # geom_text() + # print point labels
    ggtitle(title)+ggplot_theme
  
  g
}

make_multiplot <- function(variable, ylab_str, title, y_min, y_max) {
  make_subplot <- function(roi) {
    ggplot (aes(x= subj_type, y=value),data=dat[dat$variable==variable&dat$ROI==roi,])+ylab(ylab_str)+xlab('Subjects')+
      geom_boxplot(alpha=c(0.2,0.2),col=c(col1,col2),fill=c(col1,col2), lwd=0.8, width=0.4)+
      geom_point(col="gray")+ylim(y_min,y_max)+  
      ggtitle(roi)+ggplot_theme
  }
  
  a <- make_subplot('GM')
  b <- make_subplot('GL')
  c <- make_subplot('Soleus')

  g<-ggarrange(a, b, c,labels = c("a","b","c"),
             ncol = 3, nrow = 1) + ggtitle(title)
  
  g
}

make_multiplot_with_outliers <- function(variable, ylab_str, title, y_min, y_max, outlier_precision = 2, extra_ggplot_data = NULL) {
  
  make_subplot <- function(roi, extra_ggplot_data_subplot) {
    df_roi <- dat[dat$variable == variable & dat$ROI == roi, ]
    
    # Identify outliers
    df_roi$out_of_bounds <- df_roi$value < y_min | df_roi$value > y_max
    
    g <- ggplot(aes(x = subj_type, y = value), data = df_roi) + ylab(ylab_str) + xlab('Subjects') +
      geom_boxplot(alpha = c(0.2, 0.2), col = c(col1, col2), fill = c(col1, col2), lwd = 0.8, width = 0.4) +
      geom_point(col = "gray") + # ylim(y_min, y_max) +  
      ggtitle(roi) + ggplot_theme +
      scale_y_continuous(limits = c(y_min, y_max), oob = scales::squish)
    
    lim_height <- y_max - y_min
    
    g <- g + geom_segment(data = df_roi %>% filter(out_of_bounds),
                          aes(x = subj_type, y = ifelse(value < y_min, y_min+lim_height/20, y_max-lim_height/20),
                              xend = subj_type,
                              yend = ifelse(value < y_min, y_min, y_max),
                              color = subj_type),
                          arrow = arrow(type = "closed", length = unit(0.1, "inches")),
                          inherit.aes = FALSE,
                          show.legend = FALSE) +
        scale_color_manual(values = c("V" = col2, "P" = col1))
    
    # Adding text annotations for exceeded limits
    #outliers <- df_roi %>% filter(value < y_min | value > y_max)
    
    outlier_format = sprintf(' %%.%df', outlier_precision)
    
    for (s in c('V', 'P')) {
      outliers <- df_roi %>% filter(value > y_max & subj_type == s)
      outliers$index <- seq_along(outliers$value)
      
      g <- g + geom_text(data = outliers,
                         aes(x = subj_type,
                             y = ifelse(value < y_min, y_min + (index - 1) * lim_height/20, y_max - (index - 1) * lim_height/20),
                             label = sprintf(outlier_format, value),
                             hjust = ifelse(value < y_min, 1.1, -0.1)),
                         vjust = 0.5, color = rgb(0.3, 0.3, 0.3), size = 3)
      
      outliers <- df_roi %>% filter(value < y_min & subj_type == s)
      outliers$index <- seq_along(outliers$value)
      
      g <- g + geom_text(data = outliers,
                         aes(x = subj_type,
                             y = ifelse(value < y_min, y_min + (index - 1) * lim_height/20, y_max - (index - 1) * lim_height/20),
                             label = sprintf(outlier_format, value),
                             hjust = ifelse(value < y_min, 1.1, -0.1)),
                         vjust = 0.5, color = rgb(0.3, 0.3, 0.3), size = 3)
    }
    
    
    g + extra_ggplot_data_subplot
  }
  
  a <- make_subplot('GM', extra_ggplot_data[2])
  b <- make_subplot('GL', extra_ggplot_data[3])
  c <- make_subplot('Soleus', extra_ggplot_data[4])
  
  g <- ggarrange(a, b, c,
                 ncol = 3, nrow = 1) + ggtitle(title)
  
  g + extra_ggplot_data[1]
}

make_cor_plot_outliers <- function(variable1, variable2, xlims, ylims, title) {
  
  # Reshape the data
  df_wide <- dat %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    filter(subj_type %in% c("V", "P")) %>%
    filter(ROI != "Total")
  
  # identify and adjust outliers
  df_adjusted <- df_wide %>%
    mutate(
      x_adjusted = pmin(pmax(!!sym(variable1), xlims[1]), xlims[2]),
      y_adjusted = pmin(pmax(!!sym(variable2), ylims[1]), ylims[2]),
      out_of_bounds = !!sym(variable1) < xlims[1] | !!sym(variable1) > xlims[2] | 
        !!sym(variable2) < ylims[1] | !!sym(variable2) > ylims[2]
    )
  
  shape_map <- c("V" = 22, "P" = 24)
  shape_map_roi <- c("Soleus" = 22, "GM" = 23, "GL" = 24)
  color_map <- c("V" = col2, "P" = col1)
  
  g <- ggplot(df_adjusted, aes(x = x_adjusted, y = y_adjusted, shape = ROI)) +
    geom_point(aes(shape = ROI, color = subj_type), size = 3) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map_roi) +
    ggtitle(title) + ggplot_theme + xlab(variable1) + ylab(variable2)
  
  lim_width <- xlims[2] - xlims[1]
  lim_height <- ylims[2] - ylims[1]
  
  # Calculate xend and yend
  df_arrows <- df_adjusted %>% 
    filter(out_of_bounds) %>%
    mutate(
      xend = case_when(
        !!sym(variable1) < xlims[1] ~ xlims[1] - lim_width/20,
        !!sym(variable1) > xlims[2] ~ xlims[2] + lim_width/20,
        TRUE ~ x_adjusted
      ),
      yend = case_when(
        !!sym(variable2) < ylims[1] ~ ylims[1] - lim_height/20,
        !!sym(variable2) > ylims[2] ~ ylims[2] + lim_height/20,
        TRUE ~ y_adjusted
      )
    )
  
  g <- g + geom_segment(data = df_arrows,
                        aes(x = x_adjusted, y = y_adjusted,
                            xend = xend, yend = yend,
                            color = subj_type),
                        arrow = arrow(type = "closed", length = unit(0.1, "inches")),
                        inherit.aes = FALSE) 
  
  # Adding text annotations for exceeded limits
  g <- g + geom_text(
    data = df_adjusted %>% filter(!!sym(variable1) < xlims[1] | !!sym(variable1) > xlims[2]),
    aes(x = ifelse(!!sym(variable1) < xlims[1], xlims[1], xlims[2]), 
        y = y_adjusted, 
        label = sprintf('%.1f', !!sym(variable1))),
    hjust = 0, vjust = 2, color = rgb(0.3,0.3,0.3), size = 4
  ) +
    geom_text(
      data = df_adjusted %>% filter(!!sym(variable2) < ylims[1] | !!sym(variable2) > ylims[2]),
      aes(x = x_adjusted, 
          y = ifelse(!!sym(variable2) < ylims[1], ylims[1], ylims[2]), 
          label = sprintf('%.1f', !!sym(variable2))),
      hjust = 0, 
      vjust = ifelse(df_adjusted %>% filter(!!sym(variable2) < ylims[1] | !!sym(variable2) > ylims[2]) %>% pull(!!sym(variable2)) < ylims[1], 1.2, -0.2), 
      color = rgb(0.3,0.3,0.3), size = 4
    )
  
  g
}

###############################
##
## Force plot
##
###############################

g <- ggplot (aes(x=subj_type, y=force),data=force_dat)+ylab('Force [N]')+xlab('Subjects')+
  geom_boxplot(alpha=c(0.2,0.2),col=c(col1,col2),fill=c(col1,col2), lwd=0.8, width=0.3)+
  geom_point(col="gray", size=1.5)+ # geom_text() + # print point labels
  ggtitle('Maximum Voluntary Force')+ggplot_theme

save_plot(g)

###############################
##
## Plots per muscle
##
###############################


save_plot(make_multiplot_with_outliers(strain_to_use, "Strain", "Strain per ROI", 0, 0.4))
save_plot(make_multiplot_with_outliers("buildup", "Buildup rate [1/s]", "Buildup rate per ROI", 0.02, 0.1))
save_plot(make_multiplot_with_outliers("release", "Release rate [1/s]", "Release rate per ROI", 0, 0.02))
save_plot(make_multiplot_with_outliers("FF100", "Fat Fraction [%]", "FF per ROI", 0, 10, 1))
save_plot(make_multiplot_with_outliers("t2", "Water T2 [ms]", "T2 per ROI", 29.8, 40, 1))

###############################
##
## Scatter plots
##
###############################

save_plot(make_cor_plot_outliers('FF100',strain_to_use, c(0,15), c(0,0.5), 'Fat Fraction and Strain') + labs(color = 'Subject type', shape = 'Muscle') + xlab('Fat fraction [%]') + ylab('Strain'))
save_plot(make_cor_plot_outliers('t2',strain_to_use, c(25,40), c(0,0.5), 'Water T2 and Strain') + labs(color = 'Subject type', shape = 'Muscle') + xlab('Water T2 [ms]') + ylab('Strain'))
save_plot(make_cor_plot_outliers('strain_1', 'strain_2', c(0,0.5), c(-2,0), "Strain 1 and Strain 2") + labs(color = 'Subject type', shape = 'Muscle') + xlab('Strain 1') + ylab('Strain 2'))

###############################
##
## Medians and quartiles
##
###############################

odds_dataframe <- data.frame(
  Term = factor(),
  OddsRatio = numeric(),
  LowerCI = numeric(),
  UpperCI = numeric()
)

factor_levels <- c()

for(var_to_test in c(strain_to_use, 'buildup', 'release', 'FF100', 't2')) {
  for (roi in c('Soleus', 'GM', 'GL')) {
    V.quantiles = quantile(dat[dat$variable == var_to_test & dat$subj_type == 'V' & dat$ROI == roi, 'value'], c(.25,.50,.75))
    P.quantiles = quantile(dat[dat$variable == var_to_test & dat$subj_type == 'P' & dat$ROI == roi, 'value'], c(.25,.50,.75))

    df_filtered <- df_wide %>% filter(ROI == roi)
    
    # Standardize 'strain_1'
    df_filtered <- df_filtered %>%
      mutate(var_standardized = scale(!!sym(var_to_test), center = TRUE, scale = TRUE))
    
    #df_filtered <- na.omit(df_filtered)
    
    model <- glm(subj_binary ~ var_standardized, data = df_filtered, family = binomial)
    
    odds_ratios <- exp(coef(model))
    
    suppressMessages({
      # Calculate 95% confidence intervals for the model coefficients
      conf_intervals <- confint(model)
    })
    
    # Exponentiate the confidence intervals to get them for the odds ratios
    exp_conf_intervals <- exp(conf_intervals)
    
    cat(sprintf('%s\t%s\tVolunteers: %.3f (%.3f - %.3f)\tPatients: %.3f (%.3f - %.3f)\tOdds ratios: %.2f (%.2f - %.2f)\n', roi, var_to_test, V.quantiles[2], V.quantiles[1], V.quantiles[3], P.quantiles[2], P.quantiles[1], P.quantiles[3], odds_ratios[2], exp_conf_intervals[2], exp_conf_intervals[4]))
    odds_dataframe <- add_row(odds_dataframe, Term=sprintf('%s - %s', var_to_test, roi), OddsRatio = odds_ratios[2], LowerCI = exp_conf_intervals[2], UpperCI = exp_conf_intervals[4])
    factor_levels <- c(factor_levels, sprintf('%s - %s', var_to_test, roi))
  }
}

V.quantiles = quantile(force_dat[force_dat$subj_type == 'V', 'force'], c(.25,.50,.75))
P.quantiles = quantile(force_dat[force_dat$subj_type == 'P', 'force'], c(.25,.50,.75))

force_dat$subj_binary <- as.integer(force_dat$subj_type == 'P')

df_filtered <- force_dat

# Standardize 'strain_1'
df_filtered <- df_filtered %>%
  mutate(var_standardized = scale(force, center = TRUE, scale = TRUE))

#df_filtered <- na.omit(df_filtered)

model <- glm(subj_binary ~ var_standardized, data = df_filtered, family = binomial)

odds_ratios <- exp(coef(model))

suppressMessages({
  # Calculate 95% confidence intervals for the model coefficients
  conf_intervals <- confint(model)
})

# Exponentiate the confidence intervals to get them for the odds ratios
exp_conf_intervals <- exp(conf_intervals)

cat(sprintf('%s\t%s\tVolunteers: %.3f (%.3f - %.3f)\tPatients: %.3f (%.3f - %.3f)\tOdds ratios: %.2f (%.2f - %.2f)\n', '', 'Force', V.quantiles[2], V.quantiles[1], V.quantiles[3], P.quantiles[2], P.quantiles[1], P.quantiles[3], odds_ratios[2], exp_conf_intervals[2], exp_conf_intervals[4]))
odds_dataframe <- add_row(odds_dataframe, Term='Force', OddsRatio = odds_ratios[2], LowerCI = exp_conf_intervals[2], UpperCI = exp_conf_intervals[4])
factor_levels <- c(factor_levels, 'Force')

###############################
##
## Odds ratio plot
##
###############################

odds_dataframe[is.na(odds_dataframe)] = 50
odds_dataframe$Term <- factor(odds_dataframe$Term, levels = rev(factor_levels))

# Creating the plot with a logarithmic scale and a reference line at odds ratio = 1
g <- ggplot(odds_dataframe, aes(x = Term, y = OddsRatio)) +
  geom_point() +  # Plot the odds ratio as points
  geom_segment(data = odds_dataframe,
               aes(x = Term, y = LowerCI,
                   xend = Term,
                   yend = UpperCI)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # Add horizontal line at odds ratio = 1
  coord_flip() +  # Flip coordinates to make it horizontal
  theme_minimal() +  # Use a minimal theme
  ylab("Odds Ratio (Log Scale)") +  # Label for y-axis
  xlab("") +  # Remove x-axis label
  scale_y_log10(limits=c(0.05, 50)) +  # Use logarithmic scale
  ggtitle('Odds ratios')

save_plot(g)

###############################
##
## Correlations
##
###############################

for (var1 in c('FF', 't2')) {
  for(var2 in c(strain_to_use, 'buildup', 'release')) {
    correl <- cor.test(df_wide[[var1]], df_wide[[var2]], method = 'spearman', use='complete.obs')
    cat(sprintf('Correlation %s with %s: %.3f (p = %.4f)\n', var2, var1, correl$estimate, correl$p.value))
  }
}

correl <- cor.test(df_wide[['strain_1']], df_wide[['strain_2']], method = 'spearman', use='complete.obs')
cat(sprintf('Correlation strain_1 with strain_2: %.3f (p = %.4f)\n', correl$estimate, correl$p.value))