
#
# Make bar charts of orders within phyla
#
# Jason Jiang - Created: 2020/07/12
#                     Last Edit: 2021/02/12
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code for getting the orders within infected phyla, and making
# bar charts representing the orders within each phylum.
#
#
# -----------------------------------------------------------------------------


# Load useful packages
library(tidyverse)
library(readxl)


# Set working directory, this is the directory I use
setwd("P:/Shared/Microsporidia database/Figure 2 (Jason)/2C - Diversity of orders in major phyla")


# Load data containing host taxonomic information
host_taxonomy <-
  read_xlsx("Table S2 Host Masterlist.xlsx")


# We are interested in analyzing the orders within these phyla
phyla_of_interest = c('Arthropoda', 'Chordata')


# Helper functions for extracting phylum and order data from a row
get_host_phylum <- function(taxonomy_str) {
  return(strsplit(taxonomy_str, ';')[[1]][2])
}


get_host_order <- function(taxonomy_str) {
  if(lengths(strsplit(taxonomy_str, ';'))[[1]] >= 4) {
    return(strsplit(taxonomy_str, ';')[[1]][4])
  }
  else {
    return('')
  }
}


# Exclude rows w/out taxonomic data, and only consider entries containing our
# phyla of interest
# Also create new columns for the phylum and order of each species
host_taxonomy <- host_taxonomy %>%
  select(Host_formatted, Taxonomy) %>%
  filter(!is.na(Taxonomy)) %>%
  filter(grepl("Chordata", Taxonomy, fixed = TRUE) |
           grepl("Arthropoda", Taxonomy, fixed = TRUE)) %>%
  rowwise() %>%
  dplyr::mutate('Phylum' = get_host_phylum(Taxonomy),
                'Order' = get_host_order(Taxonomy)) %>%
  filter(Order != '' & Order != ' ')


# Select Arthropoda data from host_taxonomy, so we can make a frequency bar
# chart of the orders in Arthropoda
arthropoda_data <- host_taxonomy %>%
  filter(Phylum == 'Arthropoda', !grepl("Not assigned", Order, fixed = TRUE))


arthropoda_order_freqs <- as.data.frame(table(arthropoda_data$Order))


ggplot(data = arthropoda_order_freqs) +
  geom_col(mapping = aes(x = reorder(Var1, Freq), y = Freq),
           fill = 'black',
           alpha = 1, width = 0.7) +
  coord_flip() +
  labs(y = 'Frequency', x = 'Arthropoda Order') +
  theme_bw() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 24),
        axis.title.x = element_text(size = 22))


# Select Chordata data from host_taxonomy, so we can make a frequency bar chart
# of the orders in Chordata
chordata_data <- host_taxonomy %>%
  filter(Phylum == 'Chordata', !grepl("Not assigned", Order, fixed = TRUE))


chordata_order_freqs <- as.data.frame(table(chordata_data$Order))


ggplot(data = chordata_order_freqs) +
  geom_col(mapping = aes(x = reorder(Var1, Freq), y = Freq),
           fill = 'black',
           alpha = 1, width = 0.7) +
  coord_flip() +
  labs(y = 'Frequency', x = 'Chordata Order') +
  theme_bw() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 22),
        axis.title.x = element_text(size = 22))