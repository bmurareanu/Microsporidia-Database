#
# Compare spore polar tube length to spore volume.
#
# Jason Jiang - Created: 2020/10/23
#                     Last Edit: 2021/02/26
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code to analyze the relationship between spore nuclei and polar
# polar filament length.
#
#
# Thanks to Brandon for this RScript layout
#
# -----------------------------------------------------------------------------

# First, let's load in some useful packages
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("readxl")) install.packages("readxl")
if (!require("ggsignif")) install.packages("ggsignif")
if (!require("olsrr")) install.packages("olsrr")
if (!require("lmtest")) install.packages("lmtest")
if (!require("L1pack")) install.packages("L1pack")
library(tidyverse)
library(readxl)
library(ggsignif) # For adding significance bars to ggplot graphs
library(olsrr) # Tools for linear regression diagnostics
library(lmtest) # For statistical tests on regression residuals
library(L1pack) # For least absolute deviations regression
library(ggbeeswarm) # For ggplot2 compatible beeswarm plots


# Set working directory
setwd("P:/Shared/Microsporidia database/Figure 4 (Ronesh)/Figure 4C - Polar Filament and Volume")


# Open the spreadsheet as a dataframe
# This spreadsheet contains data on spore volumes, experimental polar filament
# lengths and nucleus counts.
pf_volume_and_nuclei <- read_xlsx("Fig 4E_fixed.xlsx")


# Write functions for renaming nucleus groups, and cleaning volume/
# polar filament/nucleus data for analysis
name_nucleus_group <- function(Nucleus) {
  if (is.na(Nucleus)) {
    return("No Nucleus Data")
  } else if (Nucleus == 1) {
    return("1 Nucleus")
  } else {
    return("2 Nuclei")
  }
}

clean <- function(x) {
  if(is.na(x)) {
    return(NA)
  }
  if (str_detect(x, "^\\d+?\\.")) {
    return(as.numeric(str_extract(x, "^\\d+\\.\\d+")))
  }
  return(as.numeric(str_extract(x, "^\\d+")))
}


# Great, let's reformat our data for analysis with the above functions!
pf_volume_and_nuclei <- pf_volume_and_nuclei %>%
  rowwise() %>%
  mutate(Nucleus = name_nucleus_group(clean(Nucleus)), Volume = clean(Volume),
         Experimental_PF = clean(Experimental_PF), log_Volume = log10(Volume))


# Ok, with all that data cleanup, let's start analysis.
# Let's first make box/violin plot of polar tube length across uninucleate/binucleate
# spores, to see if tube length differs across nuclei.

# FYI: this is Figure S3b
ggplot(data = filter(pf_volume_and_nuclei, Nucleus != "No Nucleus Data"),
                     aes(x = Nucleus, y = Experimental_PF,
                                       fill = Nucleus)) +
  geom_boxplot(show.legend = FALSE) +
  geom_beeswarm(show.legend = FALSE) +
  geom_signif(comparisons = list(c("1 Nucleus", "2 Nuclei")), tip_length = 0.02,
              size = 1) +
  labs(title = "Experimental PF values\nn = 68", x = "Nucleus", y = "Polar tube (um)") +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 14),
        panel.background = element_rect(fill = 'white', colour = 'black'))


# Tube lengths appears to be greater for binucleate spores.
# We'll run a wilcox/Mann-Whitney U test to see if this difference is significant.
wilcox.test(data = filter(pf_volume_and_nuclei, Nucleus != "No Nucleus Data"),
            Experimental_PF ~ Nucleus) # p = 0.017


# We'll now do linear regression for polar tube length vs volume.
# Let's first split our experimental datafrane by nucleus groups, as we've
# shown that number of nuclei is associated different polar tube lengths.
exp_1_nuc <- pf_volume_and_nuclei %>%
  filter(Nucleus == "1 Nucleus")

exp_2_nuc <- pf_volume_and_nuclei %>%
  filter(Nucleus == "2 Nuclei")


# Regression for uninucleate spores
# Use log volume for regression to help normalize residuals for regression
exp_1_nuc_ols <- lm(log_Volume ~ Experimental_PF, data = exp_1_nuc)

# Test for normal residuals
ols_plot_resid_hist(exp_1_nuc_ols)
ols_plot_resid_qq(exp_1_nuc_ols) # QQ plot looks funny, check normality with test
shapiro.test(exp_1_nuc_ols$residuals) # Fail to reject null of normal distribution

# Test for homoskedasticity
ols_plot_resid_fit(exp_1_nuc_ols)
bptest(exp_1_nuc_ols) # Fail to reject null hypothesis of constant residual variance


# Regression for binucleate spores
exp_2_nuc_ols <- lm(log_Volume ~ Experimental_PF, data = exp_2_nuc)

# Test for normal residuals
ols_plot_resid_hist(exp_2_nuc_ols)
ols_plot_resid_qq(exp_2_nuc_ols) # QQ plot seems a bit strange, so let's run statistical tests for normality
shapiro.test(exp_2_nuc_ols$residuals) # Fail to reject null

# Test for homoskedasticity
ols_plot_resid_fit(exp_2_nuc_ols)
bptest(exp_2_nuc_ols) # Fail to reject null hypothesis of constant residual variance


# Regression for all spores, regardless of nucleus data
all_ols <- lm(log_Volume ~ Experimental_PF, data = pf_volume_and_nuclei)

# Test for normal residuals
ols_plot_resid_hist(all_ols)
ols_plot_resid_qq(all_ols) # QQ plot looks funny, check normality with test
shapiro.test(all_ols$residuals) # Fail to reject null of normal distribution

# Test for homoskedasticity
ols_plot_resid_fit(all_ols)
bptest(all_ols) # Fail to reject null hypothesis of constant residual variance


# Ok, we'll do the more robust spearman correlation instead, to account for
# extreme outlying points
all_spearman_corr <- cor.test(x = pf_volume_and_nuclei$Experimental_PF,
                              y = pf_volume_and_nuclei$log_Volume,
                              method = "spearman") # p = 1.108e-12

exp_1_nuc_spearman_corr <- cor.test(x = exp_1_nuc$Experimental_PF,
                              y = exp_1_nuc$log_Volume,
                              method = "spearman") # p = 0.002248

exp_2_nuc_spearman_corr <- cor.test(x = exp_2_nuc$Experimental_PF,
                                    y = exp_2_nuc$log_Volume,
                                    method = "spearman") # p = 0.07546


################################################################################

# Perfect - let's plot the volume vs PF length scatterplots, and overlay the
# linear models.
# Make 2 plots - one for all the spores regardless of nucleus data, and one
# for spores with nucleus data only.

# Figure S3A
ggplot(data = pf_volume_and_nuclei, aes(x = Experimental_PF, y = log_Volume)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(aes(color = Nucleus)) +
  scale_color_manual(values = c("#F8766D", "#00BFC4", "black")) +
  labs(title = "All spores: R² = 0.22, p = 2.7e-09
                     \u03c1 = 0.54, p = 1.1e-12") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "top")


# Figure 4E
ggplot(data = filter(pf_volume_and_nuclei, Nucleus != "No Nucleus Data"),
       aes(x = Experimental_PF, y = log_Volume, color = Nucleus)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(title = "Uninucleate spores: R² = 0.251, p = 0.0016
                                      \u03c1 = 0.47, p = 2.2e-03
                                      
Binucleate spores: R² = 0.2226, p = 0.015
                                    \u03c1 = 0.34, p = 0.080") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "top")
