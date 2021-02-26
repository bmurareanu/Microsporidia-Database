# -----------------------------------------------------------------------------
#
# Compare spore volume to transmission
#
# Jason Jiang - Created: 2020/08/16
#                     Last Edit: 2021/02/19
#
# Reinke Lab - Microsporidia Database Project
#
# Compare spore volumes across different modes of transmission (horizontal,
# vertical, both)
#
# Thanks to Brandon for this RScript layout
#
# -----------------------------------------------------------------------------

# First, load in useful packages
library(tidyverse)
library(readxl)


# Set the working directory, this is what I use.
setwd('P:/Shared/Microsporidia database/Figure S4 (Jason)/Fig S4E - Spore volume to transmission')


# Load in the database I described earlier, of microsporidia species and their
# dimensions, volumes + infection strategies.
dimension_infection_db <- read_xlsx("NEW VOLUMES DATABASE.xlsx")


# Filter out species without transmission data
dimension_infection_db <- dimension_infection_db %>%
  filter(!is.na(Transmission) & !is.na(Volume)) %>%
  separate_rows(Volume, sep = "; ") %>%
  rowwise() %>%
  mutate(Volume_formatted = as.numeric(str_split(Volume, " ")[[1]][1])) %>%
  arrange(Volume_formatted)

# First, perform Kruskal-Wallis test to see if any median tube values are significantly
# different across microsporidia infection types
kruskal.test(data = dimension_infection_db, log10(Volume_formatted) ~ Transmission)


# Make the boxplot
ggplot(dimension_infection_db, aes(x = Transmission,
                                 y = Volume_formatted,
                                 fill = Transmission)) +
  geom_boxplot(show.legend = FALSE, lwd = 0.2) +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Transmission", y = "log_Volume", title = "\np = 0.14") +
  theme(plot.title = element_text(size = 18),
        axis.text = element_text(size = 18, color= "black"),
        axis.title = element_text(size = 18),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())