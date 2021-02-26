# -----------------------------------------------------------------------------
#
# Making bar chart of infected phyla diversity
#
# Jason Jiang - Created: 2020/07/20
#                     Last Edit: 2021/02/12
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code for creating a bar chart of the phyla infected by
# microsporidia
#
# Thanks to Brandon for this RScript layout
#
# -----------------------------------------------------------------------------


# Load useful packages
library(tidyverse)
library(readxl)


# Set working directory
setwd("P:/Shared/Microsporidia database/Figure 2 (Jason)/2A - Phyla infected by microsporidia")


# Load spreadsheet containing taxonomic data of microsporidia hosts
host_taxonomic_data <-
  read_xlsx("Table S2 Host Masterlist.xlsx")


# Narrow down dataframe to columns of interest only, and properly format all
# taxonomic data entries, and mutate a new column containing host phyla only.
host_taxonomic_data <- host_taxonomic_data %>%
  filter(!is.na(Taxonomy)) %>% # exclude hosts w/out taxonomic data
  select(Host_formatted, Taxonomy) %>%
  rowwise() %>%
  mutate(Phylum = # Mutate column with host phyla
           strsplit(Taxonomy, ';')[[1]][[2]])


# Create a dataframe of infected phyla and their frequencies.
# I think this plot will look better in Excel, so I won't be using ggplot2.
infected_phyla_freqs <-
  as.data.frame(table(host_taxonomic_data$Phylum)) %>%
  arrange(-Freq)

rownames(infected_phyla_freqs) <- c() # I don't want row names.


# Manually save this dataframe as an .xlsx file, then use Excel to create a
# frequency bar chart of infections within each phylum

# openxlsx::write.xlsx(infected_phyla_freqs,"Fig 2A - Host phyla frequencies.xlsx")