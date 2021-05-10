# -----------------------------------------------------------------------------
#
# Compare median extant species between infected and uninfected animal phyla
#
# Jason Jiang - Created: 2020/08/31
#                     Last Edit: 2021/05/03
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code for comparing the median number of extant species in infected
# and uninfected phyla, and test for significance between these medians.
#
# Thanks to Brandon for this RScript layout
#
# -----------------------------------------------------------------------------


# Load useful packages
library(plyr)
library(tidyverse)
library(readxl)
library(ggbeeswarm)
library(ggsignif)


# Set working directory
setwd("P:/Shared/Microsporidia database/Figure 2/2B - Extant species in phyla")


# Load in spreadsheet with extant species counts for infected/uninfected metazoans
extants <- read_xlsx("Phyla extant species.xlsx") %>%
  group_by(Class) %>%
  mutate(Median_sp = median(Extant_species))


# Use two-tailed Mann-Whitney U test to compare number of extant species in
# infected vs. uninfected phyla.
# Extant species in infected phyla, on average, exceed extant species in
# uninfected phyla
wilcox.test(data = extants, Extant_species ~ Class) # p = 0.000563


# Make boxplot of median extant species between infected and uninfected metazoan
# phyla.
# n = 32 metazoan phyla (15 infected, 17 uninfected)
ggplot(extants, aes(x = Class, y = Extant_species, fill = Class)) +
  scale_y_log10() +
  geom_boxplot(show.legend = FALSE) +
  geom_beeswarm(show.legend = FALSE) +
  geom_signif(size = 1.05, tip_length = 0.03, y_position=c(6.5),
              xmin=c(1), xmax=c(2), annotation=c("p = 5.6e-4"), textsize = 4.5) +
  theme_bw() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())